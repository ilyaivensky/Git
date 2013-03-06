/*
 * svm_predict.cc
 *
 *  Created on: Mar 4, 2013
 *      Author: eivensky
 */

#include "svm_predict.h"
#include "svm.h"

#include <exception>
#include <sstream>
#include <errno.h>

namespace SV {

int (*info)(const char *fmt,...) = &printf;

SVM_Predict::SVM_Predict() : 
   model(NULL),
   predict_probability_(false),  
   max_nr_attr(64),
   line(NULL)
{
}

SVM_Predict::SVM_Predict(const string & model_file, bool predict_probability) : 
    model(NULL),
    predict_probability_(false),  
    max_nr_attr(64),
    line(NULL)
{
  init(model_file, predict_probability);
}

SVM_Predict::~SVM_Predict()
{
  svm_free_and_destroy_model(&model);
}

void SVM_Predict::init(const string & model_file, bool predict_probability)
{
  predict_probability_ = predict_probability;
  
  if((model=svm_load_model(model_file.c_str()))==0)
  {
	  stringstream msg;
	  msg << "can't open model file " << model_file;
	  throw exception(msg.str().c_str());
  }
  
  model->sv_indices = NULL;
  
  if (predict_probability_)
  {
    if(svm_check_probability_model(model)==0)
      throw exception("Model does not support probability estimates"); 
  }
  else if (svm_check_probability_model(model)!=0)
    info("Model supports probability estimates, but disabled in prediction.\n");
}

double SVM_Predict::predict(const svm_node * x, double * prob_estimates) const
{
  if (predict_probability_)
  {
    int svm_type=svm_get_svm_type(model);
    if (svm_type==C_SVC || svm_type==NU_SVC)
      return svm_predict_probability(model,x,prob_estimates);
  }
   
  return svm_predict(model,x);
}

double SVM_Predict::predict(FILE *input, FILE * output)
{  
  x = (struct svm_node *) malloc(max_nr_attr*sizeof(struct svm_node));
  double retval = batch_predict(input, output);
  free(x);
  free(line);
  return retval;
}

char* SVM_Predict::readline(FILE *input)
{
  int len;

  if(fgets(line,max_line_len,input) == NULL)
    return NULL;

  while(strrchr(line,'\n') == NULL)
  {
    max_line_len *= 2;
    line = (char *) realloc(line,max_line_len);
    len = (int) strlen(line);
    if(fgets(line+len,max_line_len-len,input) == NULL)
      break;
  }
  return line;
}

double SVM_Predict::batch_predict(FILE *input, FILE *output)
{
  int correct = 0;
  int total = 0;
  double error = 0;
  double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;

  int svm_type=svm_get_svm_type(model);
  int nr_class=svm_get_nr_class(model);
  double *prob_estimates=NULL;
  int j;

  if(predict_probability_)
  {
    if (svm_type==NU_SVR || svm_type==EPSILON_SVR)
      info("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n",svm_get_svr_probability(model));
    else
    {
      int *labels=(int *) malloc(nr_class*sizeof(int));
      svm_get_labels(model,labels);
      prob_estimates = (double *) malloc(nr_class*sizeof(double));
      fprintf(output,"labels");   
      for(j=0;j<nr_class;j++)
        fprintf(output," %d",labels[j]);
      fprintf(output,"\n");
      free(labels);
    }
  }

  max_line_len = 1024;
  line = (char *)malloc(max_line_len*sizeof(char));
  while(readline(input) != NULL)
  {
    int i = 0;
    double target_label, predict_label;
    char *idx, *val, *label, *endptr;
    int inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0

    label = strtok(line," \t\n");
    if(label == NULL) // empty line
	{
		stringstream msg;
		msg << "Wrong input format at line " << total+1;
        throw exception(msg.str().c_str());
	}

    target_label = strtod(label,&endptr);
    if(endptr == label || *endptr != '\0')
	{
		stringstream msg;
		msg << "Wrong input format at line " << total+1;
        throw exception(msg.str().c_str());
	}
    
    while(1)
    {
      if(i>=max_nr_attr-1)  // need one more for index = -1
      {
        max_nr_attr *= 2;
        x = (struct svm_node *) realloc(x,max_nr_attr*sizeof(struct svm_node));
      }

      idx = strtok(NULL,":");
      val = strtok(NULL," \t");

      if(val == NULL)
        break;
      errno = 0;
      x[i].index = (int) strtol(idx,&endptr,10);
      if(endptr == idx || errno != 0 || *endptr != '\0' || x[i].index <= inst_max_index)
	  {
		  stringstream msg;
		  msg << "Wrong input format at line " << total+1;
          throw exception(msg.str().c_str());
	  }
      else
        inst_max_index = x[i].index;

      errno = 0;
      x[i].value = strtod(val,&endptr);
      if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
		   {
		  stringstream msg;
		  msg << "Wrong input format at line " << total+1;
          throw exception(msg.str().c_str());
	  }

      ++i;
    }
    x[i].index = -1;

    if (predict_probability_ && (svm_type==C_SVC || svm_type==NU_SVC))
    {
      predict_label = svm_predict_probability(model,x,prob_estimates);
      fprintf(output,"%g",predict_label);
      for(j=0;j<nr_class;j++)
        fprintf(output," %g",prob_estimates[j]);
      fprintf(output,"\n");
    }
    else
    {
      predict_label = svm_predict(model,x);
      fprintf(output,"%g\n",predict_label);
    }

    if(predict_label == target_label)
      ++correct;
    error += (predict_label-target_label)*(predict_label-target_label);
    sump += predict_label;
    sumt += target_label;
    sumpp += predict_label*predict_label;
    sumtt += target_label*target_label;
    sumpt += predict_label*target_label;
    ++total;
  }

  double retval = 0.0;

  if (svm_type==NU_SVR || svm_type==EPSILON_SVR)
  {
    info("Mean squared error = %g (regression)\n",error/total);
    info("Squared correlation coefficient = %g (regression)\n",
      ((total*sumpt-sump*sumt)*(total*sumpt-sump*sumt))/
      ((total*sumpp-sump*sump)*(total*sumtt-sumt*sumt))
      );
	retval = error/total;
  }
  else
  {
	retval = 1.0 - ((double)correct/total);
    info("Accuracy = %g%% (%d/%d) (classification)\n",
      (double)correct/total*100,correct,total);
  }
  if(predict_probability_)
    free(prob_estimates);

  return retval;
}

} // namespace

