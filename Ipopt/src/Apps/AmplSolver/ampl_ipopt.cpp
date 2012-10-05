// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "AmplTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "IpoptConfig.h"
#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif

// for printf
#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

// for parametric stuff at the end
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpMultiVectorMatrix.hpp"
#include "IpDenseVector.hpp"

using namespace Ipopt;
SmartPtr<Matrix> getSensitivityMatrix(SmartPtr<IpoptApplication> app);
SmartPtr<Vector> getDirectionalDerivative(SmartPtr<IpoptApplication> app,
					  SmartPtr<Matrix> sens_matrix);
bool doIntervalization(SmartPtr<IpoptApplication> app);
std::vector<Number> getIntervalWidths(IntervalInfoSet intervals, bool do_scaling=false);
Number Abs(Number value);


int main(int argc, char**args)
{

  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Check if executable is run only to print out options documentation
  if (argc == 2) {
    bool print_options = false;
    bool print_latex_options = false;
    if (!strcmp(args[1],"--print-options")) {
      print_options = true;
    }
    else if (!strcmp(args[1],"--print-latex-options")) {
      print_options = true;
      print_latex_options = true;
    }
    if (print_options) {
      SmartPtr<OptionsList> options = app->Options();
      options->SetStringValue("print_options_documentation", "yes");
      if (print_latex_options) {
        options->SetStringValue("print_options_latex_mode", "yes");
      }
      app->Initialize("");
      return 0;
    }
  }

  // Call Initialize the first time to create a journalist, but ignore
  // any options file
  ApplicationReturnStatus retval;
  retval = app->Initialize("");
  if (retval != Solve_Succeeded) {
    printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
    exit(-100);
  }

  // Add the suffix handler for scaling
  SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);
  // Modified for warm-start from AMPL
  suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  // Add the suffix for parameter-marking
  suffix_handler->AddAvailableSuffix("parameter", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("perturbed", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  // Add the suffix for intervall-organization
  suffix_handler->AddAvailableSuffix("intervalID", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  // Add the suffixes for branching criterion, scaling, and benefit valueing
  suffix_handler->AddAvailableSuffix("branchmode", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("scaling", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("benefit_value", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);

  SmartPtr<ParaTNLP> ampl_tnlp = new AmplTNLP(ConstPtr(app->Jnlst()),
					      app->Options(),
					      args, suffix_handler);

  // Call Initialize again to process output related options
  retval = app->Initialize();
  if (retval != Solve_Succeeded) {
    printf("ampl_ipopt.cpp: Error in second Initialize!!!!\n");
    exit(-101);
  }

  const int n_loops = 1; // make larger for profiling
  for (Index i=0; i<n_loops; i++) {
    retval = app->OptimizeTNLP(ampl_tnlp);
  }
  /*
  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<Vector> delta_s = getDirectionalDerivative(app, sens_matrix);
  if (IsValid(delta_s)) {
    delta_s->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "delta_s");
  }
*/
  doIntervalization(app);

  return 0;
}

SmartPtr<Matrix> getSensitivityMatrix(SmartPtr<IpoptApplication> app)
{
  // finalize_solution method in AmplTNLP writes the solution file

  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<const Vector> y_c = app->IpoptDataObject()->curr()->y_c();
  SmartPtr<const Vector> y_d = app->IpoptDataObject()->curr()->y_d();

  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const Matrix> opt_jac_c_p = orig_nlp->jac_c_p(*x);
  SmartPtr<const Matrix> opt_jac_d_p = orig_nlp->jac_d_p(*x);
  SmartPtr<const Matrix> opt_h_p = orig_nlp->h_p(*x, 1.0, *y_c, *y_d);
  opt_jac_c_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_c_p");
  opt_jac_d_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_d_p");
  opt_h_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_h_p");

  // Get the (factorized) KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

  // Compute ds/dp from (del K)/(del s) ds/dp = (del K)/(del p)
  Index np = orig_nlp->p()->Dim();
  SmartPtr<MultiVectorMatrixSpace> mv_space = new MultiVectorMatrixSpace(np, *x->OwnerSpace());
  SmartPtr<MultiVectorMatrix> mv = dynamic_cast<MultiVectorMatrix*>(mv_space->MakeNew());
  Number* dp_values = new Number[np];
  for (int k=0; k<np; ++k) {
    // set up current dp vector as unit vector with entry at index k
    SmartPtr<DenseVector> dp = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()))->MakeNewDenseVector();
    for (int j=0; j<np; ++j)
      dp_values[j] = 0.0;
    dp_values[k] = 1.0;
    dp->SetValues(dp_values);
    // set up iterates vector and initialize - will be rhs for linear system
    SmartPtr<IteratesVector> it_vec = app->IpoptDataObject()->curr()->MakeNewIteratesVector();
    it_vec->Set(0.0);
    // set up x part of rhs iterates vector
    SmartPtr<Vector> x_it = x->MakeNew();
    opt_h_p->MultVector(1.0, *dp, 0.0, *x_it);
    it_vec->Set_x_NonConst(*x_it);
    // set up c part of rhs iterates vector
    SmartPtr<Vector> c_it = y_c->MakeNew();
    opt_jac_c_p->MultVector(1.0, *dp, 0.0, *c_it);
    it_vec->Set_y_c_NonConst(*c_it);
    // set up d part of rhs iterates vector
    SmartPtr<Vector> d_it = y_d->MakeNew();
    opt_jac_d_p->MultVector(1.0, *dp, 0.0, *d_it);
    it_vec->Set_y_d_NonConst(*d_it);
    // do actual backsolve
    SmartPtr<IteratesVector> lhs = it_vec->MakeNewIteratesVector();
    pd_solver->Solve(1.0, 0.0, *it_vec, *lhs);
    mv->SetVector(k, *lhs->x());
  }
  delete[] dp_values;
  //  mv->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "dxdp");
  SmartPtr<Matrix> retval(GetRawPtr(mv));
  return retval;
}



SmartPtr<Vector> getDirectionalDerivative(SmartPtr<IpoptApplication> app, SmartPtr<Matrix> sens_matrix) {
  SmartPtr<const IteratesVector> curr = app->IpoptDataObject()->curr();

  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  // if perturbed values are given, compute the step and print it
  SmartPtr<const DenseVector> dp = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  SmartPtr<const DenseVectorSpace> dp_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(dp->OwnerSpace()));
  SmartPtr<DenseVector> delta_p;
  if (dp_space->HasNumericMetaData("perturbed")) {
    const std::vector<Number> perturbed = dp_space->GetNumericMetaData("perturbed");
    delta_p = dp->MakeNewDenseVector();
    const Number* dp_values = dp->Values();
    Number* new_values = new Number[dp->Dim()];
    for (int k=0; k<dp->Dim(); ++k) {
      new_values[k] = perturbed[k] - dp_values[k];
    }
    delta_p->SetValues(new_values);
    delete[] new_values;
    SmartPtr<Vector> delta_x = curr->x()->MakeNewCopy();

    sens_matrix->MultVector(1.0, *delta_p, 0.0, *delta_x);
    return delta_x;
  }
  else {
    return NULL;
  }
}


bool doIntervalization(SmartPtr<IpoptApplication> app)
{

  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<Vector> delta_s = getDirectionalDerivative(app, sens_matrix);
  if (IsValid(delta_s)) {
    //    delta_s->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "delta_s");
  }

  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));

  const Index nrows = mv_sens->NRows();
  const Index ncols = mv_sens->NCols();

  // vector to store the indexes of sensitivity matrix rows containing control related datd
  std::vector<Index> ctrl_rows;
  Number tmp_bv=0;

  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  SmartPtr<const DenseVectorSpace> p_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(p->OwnerSpace()));

  // decision variable to determine branching criterion:
  // 1 is strictly greater than, 2 is strictly smaller than, 3 is absolute value str.ly greater than
  Index branchmode = p_space->GetIntegerMetaData("branchmode")[0];
  // printf("\nbranchmode is %d \n", branchmode);

  // decision variable to enable or disable sensitivity scaling - scaling branching criterion by interval width
  bool do_scale_bc = false;
  if (p_space->GetIntegerMetaData("scaling")[0]==2)
    do_scale_bc = true;
      // if (do_scale_bc)
	// printf("\nscaling is enabled\n");
	//      else
	// printf("\nscaling is disabled\n");

  // decision variable to enable benefit value branch descision instead of simple sensitivity based decision
  // benefit value is the (scalar) product of sensitivities of both boundaries of an interval
  bool do_determine_bv = false;
  if (p_space->GetIntegerMetaData("benefit_value")[0]==2)
      do_determine_bv = true;
  //    if (do_determine_bv)
      // printf("\nbenefit value calculation is enabled\n");
  //    else
      // printf("\nbenefit value calculation is disabled\n");

  // get parameter names
  const std::vector<std::string> parnames = p_space->GetStringMetaData("idx_names");
  const Index i_p = p_space->Dim();
  std::vector<std::string> par_names_tmp;
  for (int i=0;i<i_p;i++)
    par_names_tmp.push_back(parnames[i].c_str());
  const std::vector<std::string> par_names = par_names_tmp;
  // get parameter values
  const Number* p_val = p->Values();
  std::vector<Number> par_values(i_p);
  std::copy(p_val, p_val+i_p,&par_values[0]);

  IntervalInfoSet intervals = IntervalInfoSet(p);
  intervals.printSet();

  std::vector<Number> testoutput;
  intervals.getValueVec(testoutput);
  printf("\ntestoutput.size() ist: %f\n",testoutput.size());
  for (int i=0;i<testoutput.size();i++)
    printf("\ntestoutput[%d] ist : %f\n",i,testoutput[i]);

  std::vector<Number> intervalwidths;
  intervalwidths= getIntervalWidths(intervals,do_scale_bc);

  //cycle through intervalset to assign upper and lower value vector indexes to each other
  std::vector<Index> lower_idx(int(intervals.Size()/2));
  std::vector<Index> upper_idx(int(intervals.Size()/2));
  int tmp_j =0;
  for (Index i=0;i<intervals.Size();i=i+1) {
    if (intervals.isUpper(i)) {
      intervals.getIndex(i,upper_idx[tmp_j]);
      intervals.getOtherBndIdx(upper_idx[tmp_j],lower_idx[tmp_j]);
      tmp_j++;
    }
  }

  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");

  // cycle through var space interval flags to identify and save control indexes
  for (int i=0;i<nrows;i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
    //  // printf("\ndx/dp MetaData: intervalID an der Stelle %d hat den Wert %d.\n", i, var_int_flags[i]);
  }
  std::vector<Number> benefit_values(ctrl_rows.size());
  std::vector<Index> tagged_cols(ctrl_rows.size());

  //cycle through lower/upper indexes and store the value relevant to current branching options
  for (int i =0;i<lower_idx.size();i++) {
    if (ctrl_rows.size()) {
      // get concrete sensitivity values
      std::vector<Number> s_values_lower(nrows);
      const Number* s_val = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(lower_idx[i])))->Values();
      std::copy(s_val, s_val+nrows,&s_values_lower[0]);
      std::vector<Number> s_values_upper(nrows);
      s_val = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(upper_idx[i])))->Values();
      std::copy(s_val, s_val+nrows,&s_values_upper[0]);

      for (int j=0;j<ctrl_rows.size();j++) {
	if (i==0) {
	  if (do_determine_bv) {
	    if (branchmode<3)
	      benefit_values[j]=s_values_upper[ctrl_rows.at(j)]*s_values_lower[ctrl_rows.at(j)]*intervalwidths[i];
	    else
  	      benefit_values[j]=Abs(s_values_upper[ctrl_rows.at(j)]*s_values_lower[ctrl_rows.at(j)]*intervalwidths[i]);
	    tagged_cols[j]=upper_idx[i];
	  }
	  if (branchmode<3) {
	    if (s_values_upper[ctrl_rows.at(j)]>s_values_lower[ctrl_rows.at(j)]) {
	      if (branchmode==1) {
		benefit_values[j]=s_values_upper[ctrl_rows.at(j)]*intervalwidths[i];
		tagged_cols[j]=upper_idx[i];
	      } else {
		benefit_values[j]=s_values_lower[ctrl_rows.at(j)]*intervalwidths[i];
		tagged_cols[j]=lower_idx[i];
	      }
	    } else {
	      if (branchmode==1) {
		benefit_values[j]=s_values_lower[ctrl_rows.at(j)]*intervalwidths[i];
		tagged_cols[j]=lower_idx[i];
	      } else {
		benefit_values[j]=s_values_upper[ctrl_rows.at(j)]*intervalwidths[i];
		tagged_cols[j]=upper_idx[i];
	      }
	    }
	  } else { // branchmode==3, i==0
	    if (Abs(s_values_upper[ctrl_rows.at(j)])>Abs(s_values_lower[ctrl_rows.at(j)])) {
	      benefit_values[j]=Abs(s_values_upper[ctrl_rows.at(j)]*intervalwidths[i]);
	      tagged_cols[j]=upper_idx[i];
	    } else {
	      benefit_values[j]=Abs(s_values_lower[ctrl_rows.at(j)]*intervalwidths[i]);
	      tagged_cols[j]=lower_idx[i];
	    }
	  }
	  // printf("\nbenefit_values Erstzuweisung (%e). Resultat aus: oben: %e    unten: %e   Intervalbreite: %e Spalte Obergrenze: %d Spalte Untergrenze: %d .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],upper_idx[i],lower_idx[i]);

	} else { // i=/=0
	  if (branchmode==3) {
	    if (do_determine_bv) {
	      tmp_bv = Abs(intervalwidths[i]*s_values_upper[ctrl_rows.at(j)]*s_values_lower[ctrl_rows.at(j)]);
	      // printf("\ntmp_bv is: %e   benefit_values[j] is: %e \n",tmp_bv,benefit_values[j]);
	      if (tmp_bv>benefit_values[j]) {
		benefit_values[j] = tmp_bv;
		tagged_cols[j]=upper_idx[i];
		// printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
	      }
	    } else { // no benefit_values : compare upper and lower values first, then compare the stronger one with currently stored value
	      // printf("\nupper value: %e   lower value: %e\n",s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)]);
	      if (Abs(s_values_upper[ctrl_rows.at(j)])>Abs(s_values_lower[ctrl_rows.at(j)])) {
		tmp_bv = Abs(intervalwidths[i]*s_values_upper[ctrl_rows.at(j)]);
		if (tmp_bv>benefit_values[j]) {
		  benefit_values[j] = tmp_bv;
		  tagged_cols[j]= upper_idx[i];
		  // printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
		}
	      } else { // lower value larger than upper value
		tmp_bv = Abs(intervalwidths[i]*s_values_lower[ctrl_rows.at(j)]);
		if (tmp_bv>benefit_values[j]) {
		  benefit_values[j] = tmp_bv;
		  tagged_cols[j]= lower_idx[i];
		  // printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
		}
	      }
	    }
	  } //branchmode is not 3

	  else if (do_determine_bv) {
	    tmp_bv = intervalwidths[i]*s_values_upper[ctrl_rows.at(j)]*s_values_lower[ctrl_rows.at(j)];
	    // printf("\ntmp_bv is: %e   benefit_values[j] is: %e \n",tmp_bv,benefit_values[j]);
	    if (tmp_bv>benefit_values[j]) {
	      if (branchmode ==1) {
		benefit_values[j] = tmp_bv;
		tagged_cols[j]=(upper_idx[i]);
		// printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
	      }
	    } else {
	      if (branchmode ==2) {
		benefit_values[j] = tmp_bv;
		tagged_cols[j]=(upper_idx[i]);
		// printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
	      }
	    }
	  } else { // no benefit_values : compare upper and lower values first, then compare the stronger one with currently stored value
	    // printf("\nupper value: %e   lower value: %e\n",s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)]);
	    if (s_values_upper[ctrl_rows.at(j)]>s_values_lower[ctrl_rows.at(j)]) {
	      if (branchmode==1) {
		tmp_bv = intervalwidths[i]*s_values_upper[ctrl_rows.at(j)];
	      } else { //branchmode ==2
		tmp_bv = intervalwidths[i]*s_values_lower[ctrl_rows.at(j)];
	      }
	      if (tmp_bv>benefit_values[j]) {
		if (branchmode==1){
		  benefit_values[j] = tmp_bv;
		  tagged_cols[j]= upper_idx[i];
		  // printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
		}
	      } else if (branchmode==2) {
		benefit_values[j] = tmp_bv;
		tagged_cols[j]= upper_idx[i];
		// printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
	      }
	    } else { // lower value larger than upper value
	      if (branchmode==1) {
		tmp_bv = intervalwidths[i]*s_values_lower[ctrl_rows.at(j)];
	      } else { //branchmode ==2
		tmp_bv = intervalwidths[i]*s_values_upper[ctrl_rows.at(j)];
	      }
	      if (tmp_bv>benefit_values[j]) {
		if (branchmode==1){
		  benefit_values[j] = tmp_bv;
		  tagged_cols[j]= upper_idx[i];
		  // printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
		}
	      } else if (branchmode==2) {
		benefit_values[j] = tmp_bv;
		tagged_cols[j]= upper_idx[i];
		// printf("\nbenefit_values Neuzuweisung (%e) Resultat aus: oben: %e    unten: %e   Intervalbreite: %e   tmp_bv: %e Spalte oben: %d Spalte unten: %d  .\n",benefit_values.at(j),s_values_upper[ctrl_rows.at(j)],s_values_lower[ctrl_rows.at(j)],intervalwidths[i],tmp_bv, upper_idx[i],lower_idx[i]);
	      }
	    }
	  }
	}
      }// end of j for (control cycle)
    } // end of sens size check
  }// end of sensitivity cycling for

  SmartPtr<const IteratesVector> curr = app->IpoptDataObject()->curr();
  // get critical interval and parameter for all controls
  Index tmp_idx=0;
  std::vector<Index> crit_int(tagged_cols.size());
  std::vector<Index> crit_par(tagged_cols.size());

  for (int j=0; j<tagged_cols.size();j++) {
    for (int i=0; i<i_p;i++) {
      intervals.getIndex(i,tmp_idx);
      //            // printf("\n\n aktuell untersuchter Index ist: %d\n\n",tmp_idx);
      if (tmp_idx == tagged_cols[j]) {
		// printf("\n\ncrit_int erst: %d\n\n",crit_int[j]);
	intervals.getIntervalID(i,crit_int[j]);
		// printf("\n\ncrit_int dann:: %d\n\n",crit_int[j]);
		// printf("\n\ncrit_par erst: %d\n\n",crit_par[j]);
	intervals.getParameterID(i,crit_par[j]);
		// printf("\n\ncrit_par dann: %d\n\n",crit_par[j]);
	i=i_p;
      }
    }
  }

  //write gathered information into .dat file to access with python
  std::string fname = "branch_intervals.dat";
  std::ofstream branch_intervals;
  branch_intervals.open(fname.c_str());
  char buffer[63];
  branch_intervals << "#.dat file automatically generated by AMPL intervallization routine\n#Ben Waldecker Sep 2012\n";
  for (int i=0;i<crit_int.size();i++) {
    sprintf(buffer,"\nintervalID: %d parameter: %d\n",crit_int[i],crit_par[i]);
    branch_intervals << buffer;
  }
  branch_intervals << "\n\n#end of file";
  branch_intervals.close();

  return 1;
}

IntervalInfoSet getIntInfoSet(SmartPtr<const DenseVector> parameters, const std::vector<std::string> par_names, std::vector<Number> par_values)
{
  SmartPtr<const DenseVectorSpace> p_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(parameters->OwnerSpace()));

  const std::vector<Index> intervalflags = p_space->GetIntegerMetaData("intervalID");
  const std::vector<Index> parameterflags = p_space->GetIntegerMetaData("parameter");

  //  // printf("\n\n intervalflag: %d parameterflag: %d  \n\n", intervalflags[tagged_cols[0]], parameterflags[tagged_cols[0]]);
  const Index i_p = p_space->Dim();

  // ParameterSet is to contain all parameter/interval information
  // this would prefer to be a list
  std::vector<IntervalInfo> parametersets;

  Index* tmp_par = new Index;
  Index* tmp_ID = new Index;
  bool tmp_is_upper = 0;
  IntervalInfo IntInfo;
  const std::vector<Number> p_values = par_values;

  // search for parameterentries completing one set of parameters
  for (int j =0; j< i_p; j++) {
    *tmp_par = parameterflags[j];
    *tmp_ID = intervalflags[j];
    for (int k=j+1;k<i_p;k++) {
      if (parameterflags[k] && intervalflags[k]) {
	if (*tmp_par == parameterflags[k] && *tmp_ID == intervalflags[k]) {
	  // add set to list of parametersets
	  tmp_is_upper = (par_values[j]>par_values[k]);
	  IntInfo = IntervalInfo(p_values[j],*tmp_par,*tmp_ID,j,tmp_is_upper);
	  parametersets.push_back(IntInfo);
	  IntInfo = IntervalInfo(p_values[k],*tmp_par,*tmp_ID,k,!tmp_is_upper);
	  parametersets.push_back(IntInfo);
	  k = i_p;
	}
      }
    }
  }
  IntervalInfoSet intervals = IntervalInfoSet(parametersets);
  return intervals;
}

std::vector<Number> getIntervalWidths(IntervalInfoSet intervals,bool do_scaling)
{

  std::vector<Number> par_values;
  intervals.getValueVec(par_values);

  std::vector<Number> intervalwidths(int(intervals.Size()/2));
  if (do_scaling) {

    std::vector<Index> intervalflags;
    intervals.getIntervalIDVec(intervalflags);
    // printf("\nintervalflags.size() is: %d\n",intervalflags.size());
    //    for (int i=0;i<intervalflags.size();i++)
      // printf("\nintervalflags[%d] is: %d\n",i,intervalflags[i]);

    std::vector<Index> parameterflags;
    intervals.getParameterIDVec(parameterflags);
    // printf("\nparameterflags.size() is: %d\n",intervalflags.size());
    //    for (int i=0;i<parameterflags.size();i++)
      // printf("\nparameterflags[%d] is: %d\n",i,parameterflags[i]);

    Index* tmp_par = new Index;

    // determine total parameter intervalwidths
    intervals.getParameterCount(*tmp_par);
    std::vector<Number> tmp_upper(*tmp_par);
    std::vector<Number> tmp_lower(*tmp_par);
    std::vector<bool> tmp_is_set(*tmp_par);
    for (int i=0;i<*tmp_par;i++) {
      tmp_is_set[i]=false;
    }
    for (int i=0;i<par_values.size();i++)
      printf("\npar_values[%d] ist: %f.\n",i,par_values[i]);
    for (int j=0; j<intervals.Size();j++) {
      *tmp_par = parameterflags[j]-1;
      if (!tmp_is_set[*tmp_par]) {
	tmp_upper[*tmp_par]= par_values[j];
	tmp_lower[*tmp_par]= par_values[j];
	tmp_is_set[*tmp_par]= true;
      } else {
	if (intervals.isUpper(j)) {
	  if ((tmp_upper[*tmp_par])<par_values[j]) {
	    tmp_upper[*tmp_par]=par_values[j];
	  }
	} else {
	  if ((tmp_lower[*tmp_par])>par_values[j]) {
	    tmp_lower[*tmp_par]=par_values[j];
	  }

	}
      }
    }
    std::vector<Number> total_int_widths(tmp_upper.size());
    for (int i=0;i<tmp_upper.size();i++) {
      total_int_widths[i] = tmp_upper[i]-tmp_lower[i];
    }

    //cycle through intervalset to assign upper and lower value vector indexes to each other
    std::vector<Index> lower_idx(int(intervals.Size()/2));
    std::vector<Index> upper_idx(int(intervals.Size()/2));
    int tmp_j =0;
    for (Index i=0;i<intervals.Size();i=i+1) {
      if (intervals.isUpper(i)) {
	intervals.getIndex(i,upper_idx[tmp_j]);
	intervals.getOtherBndIdx(upper_idx[tmp_j],lower_idx[tmp_j]);
	tmp_j++;
      }
    }

    // determine intervalwidths for each parametervalue - the index of the intervalwidth stored here matches the one of the parameterdata in IntervalInfoSet intervals

    for (Index i=0;i<intervalwidths.size();i++) {
      intervals.getParameterID(upper_idx[i],*tmp_par);
      intervalwidths[i]=(par_values[upper_idx.at(i)]-par_values[lower_idx.at(i)])/total_int_widths[*tmp_par-1];
      // printf("\n skalierte Intervalbreite %d für Parameter %d beträgt: %f\n",i,*tmp_par,intervalwidths[i]);
    }
  } else {
    for (int i=0;i<intervalwidths.size();i++)
      intervalwidths[i]=1;
  }
  return intervalwidths;
}

Number Abs(Number value)
{
  Number retval =0;
  if (value<0)
    retval = (-value);
  else
    retval = value;
  return retval;
}

