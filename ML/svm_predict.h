/*
 * svm_predict.h
 *
 *  Created on: Mar 4, 2013
 *      Author: eivensky
 */

#ifndef SVM_PREDICT_H_
#define SVM_PREDICT_H_

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <errno.h>
#include "libsvm-3.16\svm.h"

using namespace std;

namespace SVM {

class SVM_Predict {

public:  
  
  SVM_Predict(bool predict_probability = false);
  ~SVM_Predict();
  
  void init(const string & model_file);
  double do_predict(FILE *input, FILE * output);
  
private: 
  
  int print_null(const char *s,...) {return 0;}

  struct svm_node *x;
  int max_nr_attr;

  struct svm_model* model;
  bool predict_probability;

  char *line;
  int max_line_len;

  char* readline(FILE *input);
  
private: 
  
void exit_input_error(int line_num)
{
  fprintf(stderr,"Wrong input format at line %d\n", line_num);
  exit(1);
}

double predict(FILE *input, FILE *output);

void exit_with_help()
{
  printf(
  "Usage: svm-predict [options] test_file model_file output_file\n"
  "options:\n"
  "-b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0); for one-class SVM only 0 is supported\n"
  "-q : quiet mode (no outputs)\n"
  );
  exit(1);
}

};

} // namesapace

#endif /* SVM_PREDICT_H_ */
