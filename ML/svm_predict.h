 /*                                                                 -*- C++ -*-
 * File: svm_predict.h
 * 
 * Author: Ilya Ivensky  
 * Created on: Mar 4, 2013
 *
 * Description:
 *   C++ interface for libSVM 
 *   
 */

#ifndef SVM_PREDICT_H_
#define SVM_PREDICT_H_

#include <stdio.h>
#include <string>

using namespace std;

struct svm_node;
struct svm_model;

namespace SV {

class SVM_Predict {

public:  
  
  SVM_Predict();
  SVM_Predict(const string & model_file, 
      bool predict_probability = false);
  
  ~SVM_Predict();
  
  void init(const string & model_file, 
      bool predict_probability = false);
  
  const svm_model * get_model() const { return model; } 
  
  // Batch processing. 
  double predict(FILE *input, FILE * output);
  
  // Processing single feature set
  double predict(const svm_node * input, double * prob_estimates = NULL) const;
  
private: 
  
  double batch_predict(FILE *input, FILE *output);
  int print_null(const char *s,...) {return 0;}
  
private: 
  
  svm_node *x; // features, used only for batch processing  
  svm_model* model; // a classifier
  
  // Internal data, not passed to core
  
  bool predict_probability_;
  
  // Used to read a file input
  int max_nr_attr;
  char *line;
  int max_line_len;

  char* readline(FILE *input);
};

} // namesapace

#endif /* SVM_PREDICT_H_ */
