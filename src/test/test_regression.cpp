/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

// GSL headers.
#include <sstream>
#include <vector>
#include <cmath>

// GSL headers.
#include <gsl/gsl_matrix_double.h>

// Google Test headers.
#include "gtest/gtest.h"

// Local headers.
#include "design.hpp"
#include "gsl_fitter.hpp"
#include "regression.hpp"

//using ::testing::Eq; using ::testing::ElementsAre;
//using ::testing::DoubleNear; 
using ::testing::Test;

using std::istringstream; using std::vector;

vector<double> 
gslvec_to_stlvec(gsl_vector *gsl_vec) {
  vector<double> stl_vec;
  for (size_t ind = 0; ind < gsl_vec->size; ++ind)
    stl_vec.push_back(gsl_vector_get(gsl_vec, ind));
  return stl_vec;
}

const double max_abs_error = 0.001;

static bool doubleVectorNear(vector<double> arg, vector<double> value) {
  if (arg.size() != value.size()) 
    return false;
  
  for (size_t ind = 0; ind < arg.size(); ++ind) {
    if (fabs(arg[ind] - value[ind]) >= max_abs_error)
      return false;
  }
  
  return true;
}

static bool elementsAre(vector<size_t> arg, vector<size_t> val) {
	if (arg.size() != val.size()) return false;
	for (size_t ind = 0; ind < arg.size(); ++ind) {
		if (arg[ind] != val[ind]) return false;
	}
	return true;
}

class Unfitted2x2RegressionWithGivenParameters: public Test {
public:
  istringstream os;
  Design *design;
  Regression *regression;
  gsl_vector *parameters;
  
  void SetUp() {
    os.str("f1\tf2\ns1\t1\t1\ns2\t1\t0");
    design = new Design(os);
    regression = new Regression(*design);
    
    vector<size_t> response_total(4, 10);
		size_t t[4] = {2, 3, 7, 8};
    vector<size_t> response_meth(t, t + sizeof(t) / sizeof(t[0]));

    regression->set_response(response_total, response_meth);
    
    parameters = gsl_vector_calloc(3);
    gsl_vector_set(parameters, 0, 2.0);
    gsl_vector_set(parameters, 1, 1.5);
    gsl_vector_set(parameters, 2, 0.9);

  }
  void TearDown() {
    delete design;
    design = NULL;
    delete regression;
    regression = NULL;
    gsl_vector_free(parameters);
    parameters = NULL;
  } 
};

class Unfitted10x1Regression: public Test {
public:
  istringstream design_encoding;
  Design *design;
  Regression *regression;
  
  void SetUp() {
    design_encoding.str(
      "f\ns0\t1\ns1\t1\ns2\t1\ns3\t1\ns4\t1\ns5\t1\ns6\t1\ns7\t1\ns8\t1\ns9\t1"
      );
    design = new Design(design_encoding);
    regression = new Regression(*design);
  
    vector<size_t> response_total(10, 20);
		size_t t[10] = {0,  0,  6,  9,  7,  6,  2,  8, 10,  5};
    vector<size_t> response_meth (t, t + sizeof(t) / sizeof(t[0]));
    
    regression->set_response(response_total, response_meth);
  }

  void TearDown()  {
    delete design;
    design = NULL;
    delete regression;
    regression = NULL;
  }
};

class Fitted8x2Regression: public Test {
public:
  istringstream design_encoding;
  Design *design;
  Regression *regression;

  void SetUp() { 
    design_encoding.str("f0\tf1\ns0\t1\t0\ns1\t1\t0\ns2\t1\t0\ns3\t1\t0\n"
                                "s4\t1\t1\ns5\t1\t1\ns6\t1\t1\ns7\t1\t1");
    design = new Design(design_encoding);
    regression = new Regression(*design);
  
    vector<size_t> response_total(8, 15);
		size_t t[8] = {8, 5, 2, 2, 13, 14, 8, 5};
    vector<size_t> response_meth (t, t + sizeof(t) / sizeof(t[0]));
  
    regression->set_response(response_total, response_meth);
    
		double s[3] = {0.0, 0.0, 0.2};
    vector<double> initial_parameters (s, s + sizeof(s) / sizeof(s[0]));
    gsl_fitter(*regression, initial_parameters);
  }
  
  void TearDown() {
    delete design;
    design = NULL;
    delete regression;
    regression = NULL;
  }

};

TEST(a_regression, initializes_design_matrix) {
  istringstream os("f1\tf2\ns1\t1\t1\ns2\t1\t0");
  Design design(os);
  Regression regression(design);
	EXPECT_EQ(regression.design(), design);
}

TEST(a_regression, accepts_response) {
  istringstream os("f1\tf2\ns1\t1\t1\ns2\t1\t0");
  Design design(os);
  Regression regression(design);
  vector<size_t> response_total(4, 100);
	size_t t[4] = {64, 54, 72, 32};
  vector<size_t> response_meth (t, t + sizeof(t) / sizeof(t[0]));

  regression.set_response(response_total, response_meth);

	size_t e1 [4] = {100,100,100,100};
	size_t e2 [4] = {64,54,72,32};
	vector<size_t> ev1 (e1, e1 + sizeof(e1) / sizeof(e1[0]));
	vector<size_t> ev2 (e2, e2 + sizeof(e2) / sizeof(e2[0]));
  ASSERT_EQ(elementsAre(regression.response_total(), ev1), true);
  ASSERT_EQ(elementsAre(regression.response_meth(), ev2), true);
}

TEST_F(Unfitted2x2RegressionWithGivenParameters, computes_p) {
  EXPECT_NEAR(regression->p(0, parameters), 0.9707, max_abs_error);
}

TEST_F(Unfitted2x2RegressionWithGivenParameters, computes_loglikelihood) {
  EXPECT_NEAR(regression->loglik(parameters), -18.556, max_abs_error);
}

TEST_F(Unfitted2x2RegressionWithGivenParameters, computes_gradient) {  
  gsl_vector *grad = gsl_vector_calloc(3);
  regression->gradient(parameters, grad);
  vector<double> stl_grad = gslvec_to_stlvec(grad);
  gsl_vector_free(grad);
  
	double t[3] = {-1.77654, -0.962869, -0.935081};
  vector<double> expected (t, t + sizeof(t) / sizeof(t[0]));
  ASSERT_EQ(doubleVectorNear(stl_grad, expected), true);  
}


TEST_F(Unfitted10x1Regression, finishes_fitting_in_allotted_iterations) {
	double t[2] = {0.0, 0.1};
  vector<double> initial_parameters (t, t + sizeof(t) / sizeof(t[0])); 

	EXPECT_EQ(gsl_fitter(*regression, initial_parameters), true);
}

TEST_F(Unfitted10x1Regression, checks_number_of_initial_parameters) {
  double t[1] = {0.0}; 
  vector<double> initial_parameters (t, t + sizeof(t) / sizeof(t[0]));
  
  ASSERT_ANY_THROW(gsl_fitter(*regression, initial_parameters));
}

TEST_F(Unfitted10x1Regression, creates_initial_paramters_if_none_are_given) {
  EXPECT_EQ(gsl_fitter(*regression), true);
}

TEST(regression_10x1_design, estimates_regression_parameters_) {
  istringstream design_matrix_encoding(
    "f\ns0\t1\ns1\t1\ns2\t1\ns3\t1\ns4\t1\ns5\t1\ns6\t1\ns7\t1\ns8\t1\ns9\t1");
  Design design(design_matrix_encoding);
  Regression regression(design);
  
  vector<size_t> response_total(10, 15);
	size_t t[10] = { 0,  0,  6,  9,  7,  6,  2,  8, 10,  5 }; 
  vector<size_t> response_meth (t, t + sizeof(t) / sizeof(t[0]));
    
  regression.set_response(response_total, response_meth);
  gsl_fitter(regression);
    
  vector<double> regression_parameters = regression.fitted_parameters();
	double s[2] = {-0.679533, -1.28847}; 
  vector<double> actual (s, s + sizeof(s) / sizeof(s[0]));
  
  ASSERT_EQ(doubleVectorNear(regression_parameters, actual), true);
}


TEST(regression_8x1_design, estimates_regression_parameters) {
  istringstream design_matrix_encoding(
    "f\ns0\t1\ns1\t1\ns2\t1\ns3\t1\ns4\t1\ns5\t1\ns6\t1\ns7\t1");
  Design design(design_matrix_encoding);
  Regression regression(design);
  
  vector<size_t> response_total(8, 15);
	size_t t[8] = {2, 0, 14, 8, 15, 11, 10, 1};
  vector<size_t> response_meth (t, t + sizeof(t) / sizeof(t[0]));
  
  regression.set_response(response_total, response_meth);
  
	double s[2] = {0.0, 0.2};
  vector<double> initial_parameters (s, s + sizeof(s) / sizeof(s[0]));
  
  regression.set_response(response_total, response_meth);
  gsl_fitter(regression, initial_parameters);
    
  vector<double> regression_parameters = regression.fitted_parameters();
 
 	double u[2] = {-0.00160943, -0.02214};
  vector<double> expected (u, u + sizeof(u) / sizeof(u[0]));
  
  ASSERT_EQ(doubleVectorNear(regression_parameters, expected), true);
}

TEST_F(Fitted8x2Regression, estimates_regression_parameters) {
  vector<double> regression_parameters = regression->fitted_parameters();
	double s[3] = {-0.868078, 1.61086, -1.8563};
  vector<double> expected (s, s + sizeof(s) / sizeof(s[0]));
  ASSERT_EQ(doubleVectorNear(regression_parameters, expected), true);
}

TEST_F(Fitted8x2Regression, estimates_distribution_parameters) {      
  vector<double> distribution_parameters = 
                                  regression->fitted_distribution_parameters();
  
  EXPECT_EQ(distribution_parameters.size(), size_t(9));
  
	double s[9] = { 0.2960, 0.2960, 0.2960, 0.2960,
                  0.6776, 0.6776, 0.6776, 0.6776,
                  0.1351 };
  vector<double> expected (s, s + sizeof(s) / sizeof(s[0])); 
 
  ASSERT_EQ(doubleVectorNear(distribution_parameters, expected), true);
}

TEST_F(Fitted8x2Regression, maximum_loglikelihood) {  
  EXPECT_NEAR(regression->maximum_likelihood(), -69.8045, max_abs_error);
}

TEST_F(Fitted8x2Regression, estimates_min_methdiff) {
  EXPECT_NEAR(regression->min_methdiff(1), 0.381948, max_abs_error);
}

TEST_F(Fitted8x2Regression, estimates_log_fold_change) {
  EXPECT_NEAR(regression->log_fold_change(1), 1.61086, max_abs_error);
}
