// directly return given value if there are no different intervals
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

////////////////////////start of intervallization pseudoheader//////////////////////////////////////
class IntervalInfo : public ReferencedObject
{
public:
  IntervalInfo();
  IntervalInfo(const Number& value, const Index& parameterID, const Index& intervalID, const Index& vector_index, const bool& is_upper);
  ~IntervalInfo();
  void setParameters(const std::vector<std::string>& pnames, const std::vector<Number>& pvalues);
  void addParameter(const std::vector<std::string>& pnames, const std::vector<Number>& pvalues);
  Index getIndex() const;
  Number getValue() const;
  void setValue(const Number &value);
  Index getIntervalID() const;
  Index getParameterID() const;
  bool isUpper() const;
  void setInterval(const Index& nint);
  void printSet() const;

private:
  Number value_;
  Index parameterID_;
  Index intervalID_;
  Index index_;
  bool is_upper_;

};

class IntervalInfoSet
{
public:

  IntervalInfoSet(const std::vector<IntervalInfo>& intinfovec);
  IntervalInfoSet(SmartPtr<const DenseVector> parameters);
  ~IntervalInfoSet();
  void setIntInfoSet(const std::vector<IntervalInfo> &intinfovec);
  std::vector<IntervalInfo> getIntInfoSet() const;
  std::vector<Index> getIndexVec() const;
  Index getIndex(const Index& intindex) const;
  std::vector<Number> getValueVec() const;
  Number getValue(const Index& intindex) const;
  void setValueVec(const std::vector<Number>& valuevec);
  void setValue(const Index& intindex, const Number& value);
  std::vector<Index> getIntervalIDVec () const;
  Index getIntervalID(const Index& intindex) const;
  std::vector<Index> getIntervalIDs() const;
  std::vector<Index> getParameterIDVec() const;
  Index getParameterID(const Index& paraindex) const;
  std::vector<Index> getParameterIDs() const;
  std::vector<bool> isUpperVec() const;
  bool isUpper(const Index& isupperindex) const;
  Index getOtherBndIdx(const Index& boundindex) const;
  Index getParameterCount() const;
  Index getIntervalCount() const;
  void printSet() const;
  Index Size() const;

private:
  IntervalInfoSet();
  std::vector<IntervalInfo> intinfovec_;
  std::vector<Number> valuevec_;
  // note: at the moment, indexvec_ is obsolete, since the given vector indexes HAVE to range from 0 to # of intervalentries. if the constructor is changed towards acceptance of other indexsystems, indexvec_ turns useful again
  std::vector<Index> indexvec_;
  std::vector<Index> parameterIDvec_;
  std::vector<Index> intervalIDvec_;
  std::vector<bool> is_uppervec_;

};

class IntervalWidthScaling
{
public:
  virtual std::vector<Number> scaleIntervalWidths(const IntervalInfoSet& intervals,const Index& intervalID) const = 0;
};
IntervalWidthScaling* assignScalingMethod(SmartPtr<OptionsList> options);

class NoScaling : public IntervalWidthScaling
{
public:
  virtual std::vector<Number> scaleIntervalWidths(const IntervalInfoSet& intervals,const Index& intervalID)const;
};

class TotIntWidthScaling : public IntervalWidthScaling
{
public:
  virtual std::vector<Number> scaleIntervalWidths(const IntervalInfoSet& intervals,const Index& intervalID)const;
};

class IntWidthScaling : public IntervalWidthScaling
{
public:
  virtual std::vector<Number> scaleIntervalWidths(const IntervalInfoSet& intervals,const Index& intervalID)const;
};

struct SplitChoice
{
  Index intervalID;
  Index parameterID;
  Number reason_value;
};

struct SplitDecision
{
  Index intervalID;
  Index parameterID;
};

struct ShrinkedVector
{
  SmartPtr<Vector> vectorP;
  Index dim;
};

struct ShrinkedMatrix
{
  SmartPtr<Matrix> matrixP;
  Index ncols;
  Index nrows;
};

struct RHSVectors
{
  SmartPtr<const DenseVector> static_h_part;
  SmartPtr<const DenseVector> u_part;
  SmartPtr<const DenseVector> static_c_part;
  SmartPtr<const DenseVector> static_d_part;
  SmartPtr<const DenseVector> shift_h_part;
  SmartPtr<const DenseVector> shift_c_part;
  SmartPtr<const DenseVector> shift_d_part;
};

class BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const = 0;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const = 0;
};

BranchingCriterion* assignBranchingMethod(SmartPtr<OptionsList> options);

class RandomBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const;
};

class LargerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const;
};

class SmallerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const;
};

class AbsoluteLargerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const;
};

class Intervaluation
{
public:
  virtual SplitChoice intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<IpoptApplication> app) const= 0;
};
Intervaluation* assignIntervaluationMethod(SmartPtr<OptionsList> options);

class OneBoundIntervaluation : public Intervaluation
{
public:
  SplitChoice intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<IpoptApplication> app) const;
};

class BothBoundIntervaluation : public Intervaluation
{
public:
  SplitChoice intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<IpoptApplication> app) const;
};

class ControlSelector
{
public:
  virtual SplitDecision decideSplitControl(const std::vector<SplitChoice>& choices) const = 0;
};

class SelectFirstControl : public ControlSelector
{
public:
  SplitDecision decideSplitControl(const std::vector<SplitChoice>& choices) const;
};

ControlSelector* assignControlMethod(SmartPtr<OptionsList> options);

class ParameterShift
{
public:
  virtual RHSVectors getShiftedRHS(const Index& column,SmartPtr<IpoptApplication> app,const Index & interval) const =0;
  //  const virtual Number* getExpandedRHSValues() const =0;
  virtual SmartPtr<Matrix> getShiftedSensitivities(SmartPtr<IpoptApplication> app,const Index & interval) const =0;
  virtual RHSVectors shiftRHSVectors(const RHSVectors& originalRHS,const std::vector<Index>& intervalIDs,const Index& n_intervals, const Index& n_controls,const Index& interval) const= 0;
};

// maybe the order of deciion making should be switched - first decide for KKT lin and then for the shift
class LinearizeKKT : public ParameterShift
{
public:
  RHSVectors getShiftedRHS(const Index& column,SmartPtr<IpoptApplication> app,const Index & interval) const;
  //  const Number* getExpandedRHSValues() const;
  /* Method setting up RHS of Lin KKT as a MultiVecMatrix */
  SmartPtr<Matrix> getShiftedSensitivities(SmartPtr<IpoptApplication> app,const Index & interval) const;
  RHSVectors shiftRHSVectors(const RHSVectors& originalRHS,const std::vector<Index>& intervalIDs,const Index& n_intervals, const Index& n_controls,const Index& interval) const;
};

ParameterShift* assignShiftMethod(SmartPtr<OptionsList> options);

/////////////////////////////end of intervallization pseudo header///////////////////////////////////
SmartPtr<MultiVectorMatrix> getMVMFromMatrix(SmartPtr<Matrix> matrix);

const Index getEntryCountPerInterval(SmartPtr<const Vector>x);
Index getTotConstrCount(SmartPtr<IpoptApplication> app);
SmartPtr<Matrix> getSensitivityMatrix(SmartPtr<IpoptApplication> app);
SmartPtr<Vector> getDirectionalDerivative(SmartPtr<IpoptApplication> app,
					  SmartPtr<Matrix> sens_matrix);
bool doIntervalization(SmartPtr<IpoptApplication> app);
Number Abs(Number value);

///////////////////////////start of intervallization pseudo implementation//////////////////////////
IntervalInfo::IntervalInfo() {}

IntervalInfo::IntervalInfo(const Number& value, const Index& parameterID, const Index& intervalID, const Index& vector_index, const bool& is_upper)
{

  value_ = value;
  parameterID_ = parameterID;
  intervalID_ = intervalID;
  index_ = vector_index;
  is_upper_ = is_upper;

}

IntervalInfo:: ~IntervalInfo() {}

void IntervalInfo::setParameters(const std::vector<std::string>& pnames, const std::vector<Number>& pvalues)
{  }

void IntervalInfo::addParameter(const std::vector<std::string>& pnames, const std::vector<Number>& pvalues)
{ }

Index IntervalInfo::getIndex() const
{
  return index_;
}

Number IntervalInfo::getValue() const
{
  return value_;
}

void IntervalInfo::setValue(const Number& value)
{
  value_ = value;
}

Index IntervalInfo::getIntervalID() const
{
  return intervalID_;
}

Index IntervalInfo::getParameterID() const
{
  return parameterID_;
}

bool IntervalInfo::isUpper() const
{
  return is_upper_;
}

void IntervalInfo::setInterval(const Index& nint)
{  }

void IntervalInfo::printSet() const
{
  printf("\n value: %f parameterID: %d, intervalID: %d, index: %d, is_upper: %d \n", value_, parameterID_, intervalID_, index_, is_upper_);
}


IntervalInfoSet::IntervalInfoSet() { }

IntervalInfoSet::IntervalInfoSet(const std::vector<IntervalInfo>& intinfovec)
{
  valuevec_.clear();
  intinfovec_.clear();
  indexvec_.clear();
  parameterIDvec_.clear();
  intervalIDvec_.clear();
  is_uppervec_.clear();
  Index tmp_index;
  int i=0;

  // sort algorithm to make sure entry indexing in IntervalInfoSet matches given indexing
  while (intinfovec_.size()<intinfovec.size()) {
    tmp_index = intinfovec[i].getIndex();
    if (tmp_index==indexvec_.size()) {
      intinfovec_.push_back(intinfovec[i]);
      indexvec_.push_back(tmp_index);
      parameterIDvec_.push_back(intinfovec[i].getParameterID());
      intervalIDvec_.push_back(intinfovec[i].getIntervalID());
      is_uppervec_.push_back(intinfovec[i].isUpper());
    }
    i++;
    if (i==intinfovec.size())
      i=0;
  }

}

IntervalInfoSet::IntervalInfoSet(SmartPtr<const DenseVector> parameters)
{
  valuevec_.clear();
  intinfovec_.clear();
  indexvec_.clear();
  parameterIDvec_.clear();
  intervalIDvec_.clear();
  is_uppervec_.clear();
  SmartPtr<const DenseVectorSpace> p_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(parameters->OwnerSpace()));

  const std::vector<Index> intervalflags = p_space->GetIntegerMetaData("intervalID");
  const std::vector<Index> parameterflags = p_space->GetIntegerMetaData("parameter");
  const Index i_p = p_space->Dim();
  // get parameter values
  const Number* p_val = parameters->Values();
  std::vector<Number> par_values(i_p);
  /*std::copy(p_val, p_val+i_p,par_values);*/
  for (int i=0;i<i_p;i++)
    par_values[i] = *(p_val+i);
  valuevec_ = par_values;
  // ParameterSet is to contain all parameter/interval information
  std::vector<IntervalInfo> parametersets;

  IntervalInfo IntInfo;
  const std::vector<Number> p_values = valuevec_;
  Index* tmp_par = new Index;
  Index* tmp_ID = new Index;
  std::vector<bool> upperflags(i_p);

  // search for parameterentries completing one set of parameters
  for (int j =0; j< i_p; j++) {
    *tmp_par = parameterflags[j];
    *tmp_ID = intervalflags[j];
    for (int k=j+1;k<i_p;k++) {
      if (parameterflags[k] && intervalflags[k]) {
	if (*tmp_par == parameterflags[k] && *tmp_ID == intervalflags[k]) {
	  upperflags[j] = (valuevec_[j]>valuevec_[k]);
	  upperflags[k] = (!upperflags[j]);
	  k = i_p;
	}
      }
    }
    IntInfo = IntervalInfo(p_values[j],*tmp_par,*tmp_ID,j,upperflags[j]);
    intinfovec_.push_back(IntInfo);
    indexvec_.push_back(j);
    parameterIDvec_.push_back(*tmp_par);
    intervalIDvec_.push_back(*tmp_ID);
    is_uppervec_.push_back(upperflags[j]);
  }

}

IntervalInfoSet::~IntervalInfoSet() {}

void IntervalInfoSet::setIntInfoSet(const std::vector<IntervalInfo> &intinfovec)
{
  intinfovec_=intinfovec;
}

std::vector<IntervalInfo> IntervalInfoSet::getIntInfoSet() const
{
  return intinfovec_;
}

std::vector<Index> IntervalInfoSet::getIndexVec() const
{
  return indexvec_;
}

Index IntervalInfoSet::getIndex(const Index& intindex) const
{
  if (intindex<indexvec_.size())
    return indexvec_[intindex];
  else {
    return -1;
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::getIndex() call with out of range index!\n");
  }
}

std::vector<Number> IntervalInfoSet::getValueVec() const
{
  return valuevec_;
}

Number IntervalInfoSet::getValue(const Index& intindex) const
{
 return intinfovec_[intindex].getValue();
}

void IntervalInfoSet::setValueVec(const std::vector<Number>& valuevec)
{
  valuevec_ = valuevec;
}

void IntervalInfoSet::setValue(const Index& intindex,const Number& value)
{
  intinfovec_[intindex].setValue(value);
}

std::vector<Index> IntervalInfoSet::getIntervalIDVec() const
{
  return intervalIDvec_;
}

Index IntervalInfoSet::getIntervalID(const Index& intindex) const
{
  if (intindex<intervalIDvec_.size())
    return intervalIDvec_[intindex];
  else {
    return -1;
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::getIntervalID() call with out of range index!\n");
  }
}

std::vector<Index> IntervalInfoSet::getIntervalIDs() const
{
  std::vector<Index> tmp_IDs;
  Index nints = this->getIntervalCount();
  Index* tmp_new_ID = new Index;
  Index* tmp_ID = new Index;
  *tmp_new_ID = intervalIDvec_[0];
  // assign first value
  for (int i=0;i<intervalIDvec_.size();i++)
    if (intervalIDvec_[i]<*tmp_new_ID)
      *tmp_new_ID=intervalIDvec_[i];
  tmp_IDs.push_back(*tmp_new_ID);

  //assign follow up values
  while (tmp_IDs.size()<nints) {
    *tmp_ID = *tmp_new_ID;
    for (int i=0;i<intervalIDvec_.size();i++) {
      if (intervalIDvec_[i]>*tmp_ID && (*tmp_ID==*tmp_new_ID))
 	*tmp_new_ID=intervalIDvec_[i];
      if (intervalIDvec_[i]>*tmp_ID && intervalIDvec_[i]<*tmp_new_ID)
	*tmp_new_ID = intervalIDvec_[i];
    }
    tmp_IDs.push_back(*tmp_new_ID);
  }
  return tmp_IDs;
}

std::vector<Index> IntervalInfoSet::getParameterIDVec() const
{
  return parameterIDvec_;
}

Index IntervalInfoSet::getParameterID(const Index& paraindex) const
{
  if (paraindex<parameterIDvec_.size())
    return parameterIDvec_[paraindex];
  else {
    return -1;
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::getParameterID() call with out of range index!\n");
  }
}

std::vector<Index> IntervalInfoSet::getParameterIDs() const
{
  std::vector<Index> tmp_IDs;
  Index npars = this->getParameterCount();
  Index* tmp_new_ID = new Index;
  Index* tmp_ID = new Index;
  *tmp_new_ID = parameterIDvec_[0];
  // assign first value
  for (int i=0;i<parameterIDvec_.size();i++)
    if (parameterIDvec_[i]<*tmp_new_ID)
      *tmp_new_ID=parameterIDvec_[i];
  tmp_IDs.push_back(*tmp_new_ID);

  //assign follow up values
  while (tmp_IDs.size()<npars) {
    *tmp_ID = *tmp_new_ID;
    for (int i=0;i<parameterIDvec_.size();i++) {
      if (parameterIDvec_[i]>*tmp_ID && (*tmp_ID==*tmp_new_ID))
 	*tmp_new_ID=parameterIDvec_[i];
      if (parameterIDvec_[i]>*tmp_ID && parameterIDvec_[i]<*tmp_new_ID)
	*tmp_new_ID = parameterIDvec_[i];
    }
    tmp_IDs.push_back(*tmp_new_ID);
  }
  return tmp_IDs;
}

std::vector<bool> IntervalInfoSet::isUpperVec() const
{
  return is_uppervec_;
}

bool IntervalInfoSet::isUpper(const Index& isupperindex) const
{
  if (isupperindex<is_uppervec_.size())
    return is_uppervec_[isupperindex];
  else {
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::isUpper() call with out of range index!\n");
    printf("\nsize of is_uppervec_: %d \n", is_uppervec_.size());
    return 0;
  }
}

Index IntervalInfoSet::getOtherBndIdx(const Index& boundindex) const
{
  for (int i=0;i<intinfovec_.size();i++){
    if (parameterIDvec_[boundindex]==parameterIDvec_[i] && intervalIDvec_[boundindex]==intervalIDvec_[i] && is_uppervec_[boundindex]!=is_uppervec_[i]){
      return indexvec_[i];
    }
  }
  return 0;
}

Index IntervalInfoSet::getParameterCount() const
{
  Index tmp_count =0;
  for (int i=0;i<parameterIDvec_.size();i++) {
    if (i==0)
      tmp_count = parameterIDvec_[i];
    if (parameterIDvec_[i]>tmp_count)
      tmp_count = parameterIDvec_[i];
  }
  return tmp_count;
}

Index IntervalInfoSet::getIntervalCount() const
{
  Index tmp_count= 0;
  for (int i=0;i<intervalIDvec_.size();i++) {
    if (i==0)
      tmp_count = intervalIDvec_[i];
    if (intervalIDvec_[i]>tmp_count)
      tmp_count = intervalIDvec_[i];
  }
  return tmp_count;
}

void IntervalInfoSet::printSet() const
{
  for (int i=0; i<intinfovec_.size();i++){
    printf("\n\nIntervalInfoSet Eintrag %d:\n", i);
    intinfovec_[i].printSet();
    printf("\n");
  }
}

Index IntervalInfoSet::Size() const
{
  return intinfovec_.size();
}

IntervalWidthScaling* assignScalingMethod(SmartPtr<OptionsList> options)
{
  std::string scalingmode;
  if (options->GetStringValue("scalingmode",scalingmode ,"")){
    if (scalingmode =="none")
      return new NoScaling();
     else if (scalingmode=="total_interval_widths")
       return new TotIntWidthScaling();
     else if (scalingmode=="interval_widths")
       return new IntWidthScaling();
  } else {
    printf("\nNo scalingmode given!\n");
  }
    return new NoScaling();
}


std::vector<Number> NoScaling::scaleIntervalWidths(const IntervalInfoSet& intervals,const Index& intervalID) const
{
  // no scaling: just return a vector with apropriate length containing 1s
  Index npars = intervals.getParameterCount();
  std::vector<Number> retval(npars);
  for (int i=0;i<npars;i++)
    retval[i]=1;
  return retval;
}

std::vector<Number> TotIntWidthScaling::scaleIntervalWidths(const IntervalInfoSet& intervals,const Index& intervalID) const
{
  Index npars = intervals.getParameterCount();
  std::vector<Number> retval(npars);
  for (int i=0;i<npars;i++)
    retval[i]=0;
  std::vector<Number> par_values = intervals.getValueVec();
  std::vector<Index> intervalflags = intervals.getIntervalIDVec();
  std::vector<Index> parameterflags = intervals.getParameterIDVec();
  std::vector<Index> parameterIDs = intervals.getParameterIDs();
  //store largest and smallest interval bounds for each parameter
  std::vector<Number> tmp_largest(npars);
  std::vector<Number> tmp_smallest(npars);
  std::vector<bool> tmp_is_set(npars);
  for (int i=0;i<npars;i++) {
    tmp_is_set[i] = false;
  }
  // determine (total) parameter intervalwidths
  for (int i=0;i<npars;i++) {
    for (int j=0;j<intervals.Size();j++) {
      if (parameterflags[j] == parameterIDs[i]) {
	// largest and smallest values are interval independent
	if (!tmp_is_set[i]) {
	  tmp_largest[i] = par_values[j];
	  tmp_smallest[i] = par_values[j];
	  tmp_is_set[i]= true;
	}
	else {
	  if (par_values[j] < tmp_smallest[i])
	    tmp_smallest[i] = par_values[j];
	  if (par_values[j] > tmp_largest[i])
	    tmp_largest[i] = par_values[j];
	}
	// real interval width for each parameter only needs to be determined for the intervalID in question
	if (intervalflags[j] == intervalID) {
	  if (intervals.isUpper(j))
	    retval[i] += par_values[j];
	  else
	    retval[i] -= par_values[j];
	}
      }
    }
  }
  // total interval width for each parameter is difference of largest and smallest values
  std::vector<Number> total_int_widths(npars);
  for (int i=0;i<npars;i++) {
    total_int_widths[i] = tmp_largest[i]-tmp_smallest[i];
    retval[i] = retval[i] / total_int_widths[i];
  }
  return retval;
}

std::vector<Number> IntWidthScaling::scaleIntervalWidths(const IntervalInfoSet& intervals,const Index& intervalID) const
{
  Index npars = intervals.getParameterCount();
  std::vector<Number> retval(npars);
  for (int i=0;i<npars;i++)
    retval[i]=0;
  std::vector<Number> par_values = intervals.getValueVec();
  std::vector<Index> intervalflags = intervals.getIntervalIDVec();
  std::vector<Index> parameterflags = intervals.getParameterIDVec();
  std::vector<Index> parameterIDs = intervals.getParameterIDs();

  for (int i=0;i<npars;i++) {
    for (int j=0;j<intervals.Size();j++) {
      if (parameterflags[j] == parameterIDs[i] & intervalflags[j] == intervalID) {
	// real interval width for each parameter only needs to be determined for the intervalID in question
	if (intervals.isUpper(j))
	  retval[i] += par_values[j];
	else
	  retval[i] -= par_values[j];
      }
    }
  }

  return retval;
}

BranchingCriterion* assignBranchingMethod(SmartPtr<OptionsList> options)
{
  std::string branchmode;
  if (options->GetStringValue("branchmode",branchmode ,"")){
    if (branchmode =="smallest")
      return new SmallerBranching();
     else if (branchmode=="largest")
       return new LargerBranching();
     else if (branchmode=="random")
       return new RandomBranching();
     else if (branchmode=="abs_largest")
       return new AbsoluteLargerBranching();
  } else {
    printf("\nNo branchmode given!\n");
  }
    return new RandomBranching();
}

std::vector<SplitChoice> RandomBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const
{
  // RandomBranching not implemented yet.
}

SplitChoice RandomBranching::branchInterval(const std::vector<SplitChoice>& splitchoices) const
{
  // RandomBranching not implemented yet.
}

std::vector<SplitChoice>  LargerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const
{

  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  IntervalInfoSet intervals = IntervalInfoSet(p);
  SmartPtr<OptionsList> options = app->Options();

// cycle through var space interval flags to identify and save control indexes
  const Index nrows = mv_sens->NRows();
  std::vector<Index> ctrl_rows;
  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");
  for (int i=0;i<nrows;i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
  }
  printf("\nLargerBranching::branchSensitivityMatrix: ctrl_rows.size() is: %d\n",ctrl_rows.size());
  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs = intervals.getIntervalIDs();
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater = assignIntervaluationMethod(options);
  assert(intervalIDs.size());
  for (Index i=0;i<ctrl_rows.size();i++) {
    //with only one control and only one interval, the only decision left is which parameter to split. thats done by intervaluation
    retval.push_back(intervaluater->intervaluateInterval(ctrl_rows[i],intervalIDs[0],app));
    if (intervalIDs.size()>1) {
      //with more than one interval, the branchingalgorithm is to pick the critical interval and return its data
      for (Index j=1;j<intervalIDs.size();j++) {
	tmp_split_choice = intervaluater->intervaluateInterval(ctrl_rows[i],intervalIDs[j],app);
	if (tmp_split_choice.reason_value>retval[i].reason_value)
	  retval[i] = tmp_split_choice;
      }
    }
  }
  return retval;
}

SplitChoice LargerBranching::branchInterval(const std::vector<SplitChoice>& splitchoices) const
{
  assert(splitchoices.size());
  // directly return given value if there are no different intervals
  SplitChoice retval=splitchoices[0];
  // otherwise pick the interval critical to branchmode
  if (splitchoices.size()>1) {
    for (int i=1;i<splitchoices.size();i++) {
      if (splitchoices[i].reason_value>retval.reason_value)
	retval=splitchoices[i];
    }
  }
  return retval;
}

std::vector<SplitChoice> SmallerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const
{

  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  IntervalInfoSet intervals = IntervalInfoSet(p);
  SmartPtr<OptionsList> options = app->Options();

// cycle through var space interval flags to identify and save control indexes
  const Index nrows = mv_sens->NRows();
  std::vector<Index> ctrl_rows;
  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");
  for (int i=0;i<nrows;i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
  }
  printf("\nSmallerBranching::branchSensitivityMatrix: ctrl_rows.size() is: %d\n",ctrl_rows.size());
  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs = intervals.getIntervalIDs();
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater = assignIntervaluationMethod(options);
  assert(intervalIDs.size());
  for (Index i=0;i<ctrl_rows.size();i++) {
    assert (intervalIDs.size());
    //with only one control and only one interval, the only decision left is which parameter to split. thats done by intervaluation
    retval.push_back(intervaluater->intervaluateInterval(ctrl_rows[i],intervalIDs[0],app));
    if (intervalIDs.size()>1) {
      //with more than one interval, the branchingalgorithm is to pick the critical interval and return its data
      for (Index j=1;j<intervalIDs.size();j++) {
	tmp_split_choice = intervaluater->intervaluateInterval(ctrl_rows[i],intervalIDs[j],app);
	if (tmp_split_choice.reason_value<retval[i].reason_value)
	  retval[i] = tmp_split_choice;
      }
    }
  }
  return retval;
}

SplitChoice SmallerBranching::branchInterval(const std::vector<SplitChoice>& splitchoices) const
{
  assert(splitchoices.size());
  // directly return given value if there are no different intervals
  SplitChoice retval=splitchoices[0];
  // otherwise pick the interval critical to branchmode
  if (splitchoices.size()>1) {
    for (int i=1;i<splitchoices.size();i++) {
      if (splitchoices[i].reason_value<retval.reason_value) {
	retval=splitchoices[i];
      }
    }
  }
  return retval;
}

std::vector<SplitChoice> AbsoluteLargerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const
{

  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  IntervalInfoSet intervals = IntervalInfoSet(p);
  SmartPtr<OptionsList> options = app->Options();

// cycle through var space interval flags to identify and save control indexes
  const Index nrows = mv_sens->NRows();
  std::vector<Index> ctrl_rows;
  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");
  for (int i=0;i<nrows;i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
  }

  printf("\nAbsoluteLargerBranching::branchSensitivityMatrix: ctrl_rows.size() is: %d\n",ctrl_rows.size());

  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs = intervals.getIntervalIDs();
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater = assignIntervaluationMethod(options);
  assert(intervalIDs.size());
  for (Index i=0;i<ctrl_rows.size();i++) {
    //with only one control and only one interval, the only decision left is which parameter to split. thats done by intervaluation
    retval.push_back(intervaluater->intervaluateInterval(ctrl_rows[i],intervalIDs[0],app));
    if (intervalIDs.size()>1) {
      //with more than one interval, the branchingalgorithm is to pick the critical interval and return its data
      for (Index j=1;j<intervalIDs.size();j++) {
	tmp_split_choice = intervaluater->intervaluateInterval(ctrl_rows[i],intervalIDs[j],app);
	if (Abs(tmp_split_choice.reason_value)>Abs(retval[i].reason_value))
	  retval[i] = tmp_split_choice;
      }
    }
  }
  return retval;
}

SplitChoice AbsoluteLargerBranching::branchInterval(const std::vector<SplitChoice>& splitchoices) const
{
  assert(splitchoices.size());
  // directly return given value if there are no different intervals
  SplitChoice retval=splitchoices[0];
  // otherwise pick the interval critical to branchmode
  if (splitchoices.size()>1) {
    for (int i=1;i<splitchoices.size();i++) {
      if (splitchoices[i].reason_value>retval.reason_value)
	retval=splitchoices[i];
    }
  }
  return retval;
}

Intervaluation* assignIntervaluationMethod(SmartPtr<OptionsList> options)
{
  std::string branchvalue;
  if (options->GetStringValue("branchvalue",branchvalue ,"")){
    if (branchvalue =="bound")
      return new OneBoundIntervaluation();
     else if (branchvalue=="product")
       return new BothBoundIntervaluation();
  } else {
    printf("\nNo branchvalue given!\n");
  }
    return new OneBoundIntervaluation();

}

SplitChoice OneBoundIntervaluation::intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<IpoptApplication> app) const
{
  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  IntervalInfoSet intervals = IntervalInfoSet(p);
  SmartPtr<OptionsList> options = app->Options();

  SplitChoice retval;
  BranchingCriterion* int_brancher = assignBranchingMethod(options);
  IntervalWidthScaling* scalingmode = assignScalingMethod(options);
  std::vector<Number> intervalwidths =  scalingmode->scaleIntervalWidths(intervals,intervalID);
  std::vector<Number> p_values = intervals.getValueVec();
  std::vector<Index> ID_vec = intervals.getIntervalIDVec();
  std::vector<Index> para_vec = intervals.getParameterIDVec();
  Index npars = intervals.getParameterCount();
  std::vector<SplitChoice> int_splitchoices(2*npars);
  std::vector<Index> parameterIDs = intervals.getParameterIDs();
  Number tmp_value;
  Index tmp_idx;
  //cycle through the different parameter values for the given interval and determine branching objective values
  for (int i=0; i<npars;i++) {
    int_splitchoices[2*i].parameterID = parameterIDs[i];
    int_splitchoices[2*i].intervalID = intervalID;
    int_splitchoices[2*i+1].parameterID = parameterIDs[i];
    int_splitchoices[2*i+1].intervalID = intervalID;
    for (int j=0;j<intervals.Size();j++) {
      if (ID_vec[j]==intervalID && para_vec[j]==parameterIDs[i]) {
      // parameter entry for interval in question found!
	tmp_idx = intervals.getOtherBndIdx(j);
	const Number* sensitivity_value_1 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(j)))->Values();
	const Number* sensitivity_value_2 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(tmp_idx)))->Values();
	int_splitchoices[2*i].reason_value = *(sensitivity_value_1+controlrow)*intervalwidths[i];
	int_splitchoices[2*i+1].reason_value = *(sensitivity_value_2+controlrow)*intervalwidths[i];
	j=intervals.Size();
      }
    }
  }
  //let intbrancher decide, which parameter is critical on this interval
  retval = int_brancher->branchInterval(int_splitchoices);
  return retval;
}

SplitChoice BothBoundIntervaluation::intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<IpoptApplication> app) const
{
  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  IntervalInfoSet intervals = IntervalInfoSet(p);
  SmartPtr<OptionsList> options = app->Options();

  SplitChoice retval;
  BranchingCriterion* int_brancher = assignBranchingMethod(options);
  IntervalWidthScaling* scalingmode = assignScalingMethod(options);
  std::vector<Number> intervalwidths =  scalingmode->scaleIntervalWidths(intervals,intervalID);
  std::vector<Number> p_values = intervals.getValueVec();
  std::vector<Index> ID_vec = intervals.getIntervalIDVec();
  std::vector<Index> para_vec = intervals.getParameterIDVec();
  Index npars = intervals.getParameterCount();
  std::vector<SplitChoice> int_splitchoices(2*npars);
  std::vector<Index> parameterIDs = intervals.getParameterIDs();
  Number tmp_value;
  Index tmp_idx;
  //cycle through the different parameter values for the given interval and determine branching objective values
  for (int i=0; i<npars;i++) {
    int_splitchoices[i].parameterID = parameterIDs[i];
    int_splitchoices[i].intervalID = intervalID;
    for (int j=0;j<intervals.Size();j++) {
      if (ID_vec[j]==intervalID && para_vec[j]==parameterIDs[i]) {
      // parameter entry for interval in question found!
	tmp_idx = intervals.getOtherBndIdx(j);
	const Number* sensitivity_value_1 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(j)))->Values();
	const Number* sensitivity_value_2 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(tmp_idx)))->Values();
	tmp_value = *(sensitivity_value_1+controlrow);
	int_splitchoices[i].reason_value = tmp_value*(*(sensitivity_value_2+controlrow))*intervalwidths[i];
	j=intervals.Size();
      }
    }
  }
  //let intbrancher decide, which parameter is critical on this interval
  retval = int_brancher->branchInterval(int_splitchoices);
  return retval;
}

SplitDecision SelectFirstControl::decideSplitControl(const std::vector<SplitChoice>& choices) const
{
  assert(choices.size());
  //simply returns first entry of SplitDecisions
  SplitDecision retval;
  retval.intervalID = choices[0].intervalID;
  retval.parameterID = choices[0].parameterID;
  return retval;
}

ControlSelector* assignControlMethod(SmartPtr<OptionsList> options)
{
  // no other Control methods implemented yet!
  return new SelectFirstControl();
}

const Index getEntryCountPerInterval(SmartPtr<const Vector> x)
{
  SmartPtr<const DenseVectorSpace> x_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x->OwnerSpace()));
  //  printf("\ngetEntryCountPerInterval: Vector Dimension is %d\n",x_space->Dim());
  if (!x_space->HasIntegerMetaData("intervalID"))
    return 0;
  const std::vector<Index> intervalflags =  x_space->GetIntegerMetaData("intervalID");
  const std::vector<std::string> entrynames =  x_space->GetStringMetaData("idx_names");
  assert(intervalflags.size());
  assert(entrynames.size());
  Index tmp_ID =0;
  Index retval = 0;
  const Index n_x = x_space->Dim();
  // assign first non-zero value to tmp_ID
  for (Index i=0;i<n_x;i++) {
    if (!tmp_ID && intervalflags[i]) {
      tmp_ID=intervalflags[i];
      retval++;
    } else {
      if (tmp_ID == intervalflags[i]) {
	retval++;
      }
    }
  }
  //    printf("\ngetEntryCountPerInterval: EntryCount is is %d\n",retval);
  return retval;
}


SmartPtr<Vector> getVectorFromMatrix(const Index& colindex, SmartPtr<Matrix>& matrix)
{

}

SmartPtr<Vector> extractVecFromMat(const Index& colindex, SmartPtr<Matrix>& matrix, const std::vector<Index>& indexlist)
{/*
  // get vector from which to extract target values
  SmartPtr<DenseVector> retvec = dynamic_cast<DenseVector*>(GetRawPtr(getVectorFromMatrix(colindex,matrix)));
  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(indexlist.size());
  SmartPtr<DenseVector> retval = retval_space->MakeNewDenseVector();
  const Number* rr_values = retvec->Values();
  // check if original vector is long enough for the operation:
  Index ori_length = (retvec->OwnerSpace())->Dim();
  assert(ori_length>indexlist.size());

  std::vector<const Number*> tmp_values(indexlist.size());
   for (Index i=0;i<indexlist.size();i++) {
     tmp_values[i]=(rr_values+i);
   }
  retval->SetValues(rr_values);
  // mal checken, ob der Umweg über tmp_... weggelassen werden kann!
 */
}

SmartPtr<MultiVectorMatrix> getMVMFromMatrix(SmartPtr<const Matrix> matrix)
{
  // set up neccessary Spaces and return value candidates
  const Index n_cols = (matrix->OwnerSpace())->NCols();
  const Index n_rows = (matrix->OwnerSpace())->NRows();
  //  printf("\nampl_ipopt: getMVMFromMatrix: original matrix size: NROWS: %d, NCOLS: %d.",n_rows,n_cols);
  SmartPtr<DenseVectorSpace> retvec_space = new DenseVectorSpace(n_rows);
  SmartPtr<DenseVectorSpace> ret_space = new DenseVectorSpace(n_cols);
  SmartPtr<MultiVectorMatrixSpace> retval_space = new MultiVectorMatrixSpace(n_cols,*retvec_space);
  SmartPtr<MultiVectorMatrix> retval = retval_space->MakeNewMultiVectorMatrix();
  retval->FillWithNewVectors();

  Number* unit_values = new Number[n_cols];
  for (int i=0;i<n_cols;i++) {
    SmartPtr<DenseVector> retvec = dynamic_cast<const DenseVector*>(retvec_space->MakeNewDenseVector());
    SmartPtr<DenseVector> unitvector = dynamic_cast<const DenseVector*>(ret_space->MakeNewDenseVector());
    for (int j=0;j<n_cols;j++)
      unit_values[j]=0.0;
    unit_values[i] = 1.0;
    unitvector->SetValues(unit_values);
    matrix->MultVector(1.0,*unitvector,0.0,*retvec);
    retval->SetVector(i,*retvec);
  }
  return retval;
}

  //const Number* LinearizeKKT::getExpandedRHSValues(SmartPtr<DenseVector>h,SmartPtr<DenseVector>c,SmartPtr<DenseVector>d,const Index& intervalID, const std::vector<Index>& intervalIDs,const Index& n_p, const Index& n_i) const
/*
{
  const std::vector<Index> intervalIDs = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_space))->GetIntegerMetaData("intervalID");
  Index curr_pos =0;
  SmartPtr<const DenseVector> h;
  SmartPtr<const DenseVector> c_ori_vec;
  SmartPtr<const DenseVector> d_ori_vec;
  Number* retvalues = new Number[n_rows];
  Index* skip_h_indices = new Index[xperinterval];
  Index* skip_c_indices = new Index[cperinterval];
  Index* skip_d_indices = new Index[dperinterval];
  Index* skip_u_indices = new Index[n_u];
  Index skip_cnt = 0;
  Index u_cnt = 0;

  // fill rhs_matrix with values
  for (int i=0;i<n_p;i++) {
    SmartPtr<DenseVector> retvec = dynamic_cast<const DenseVector*>(retvec_space->MakeNewDenseVector());
    h = dynamic_cast<const DenseVector*>(GetRawPtr(h_MVM->GetVector(i)));
    const Number* h_values = h->Values();
    c = dynamic_cast<const DenseVector*>(GetRawPtr(c_MVM->GetVector(i)));
    const Number* c_values = c->Values();
    d = dynamic_cast<const DenseVector*>(GetRawPtr(d_MVM->GetVector(i)));
    const Number* d_values = d->Values();

    // add objective related data
    for (int j=0;j<x->Dim();j++) {
      if (intervalIDs[j]==interval) {
	skip_h_indices[skip_cnt]=j;
	skip_cnt++;
      } else if (!intervalIDs[j]) {
	skip_u_indices[u_cnt]=j;
	u_cnt++;
      } else {
	retvalues[curr_pos]=h_ori_values[j];
	curr_pos++;
      }
    }
    skip_cnt=0;
    // add equality related data
    for (int j=0;j<rhs_c->NRows();j++) {
      if ((interval-1)*cperinterval<=j && j<interval*cperinterval) {
	skip_c_indices[skip_cnt]=j;
	skip_cnt++;
      }	else {
	retvalues[curr_pos]=c_ori_values[j];
	curr_pos++;
      }
    }
    skip_cnt=0;
    // add inequality related data
    for (int j=0;j<rhs_d->NRows();j++) {
      if ((interval-1)*dperinterval<=j && j<interval*dperinterval) {
	skip_d_indices[skip_cnt]=j;
	skip_cnt++;
      }	else {
	retvalues[curr_pos]=d_ori_values[j];
	curr_pos++;
      }
    }

    // add the duplicated interval values (objective, equalities, inequalities)
    for (int j=0;j<xperinterval;j++) {
      retvalues[curr_pos]= h_ori_values[skip_h_indices[j]];
      curr_pos++;
    }
    for (int j=0;j<cperinterval;j++) {
      retvalues[curr_pos]= c_ori_values[skip_c_indices[j]];
      curr_pos++;
    }
    for (int j=0;j<dperinterval;j++) {
      retvalues[curr_pos]= d_ori_values[skip_d_indices[j]];
      curr_pos++;
    }

    // add control related data
    for (int j=0;j<n_u;j++) {
      retvalues[curr_pos]= h_ori_values[skip_u_indices[j]];
      curr_pos++;
    }

    // add the duplicated interval values (objective, equalities, inequalities)
    for (int j=0;j<xperinterval;j++) {
      retvalues[curr_pos]= h_ori_values[skip_h_indices[j]];
      curr_pos++;
    }
    for (int j=0;j<cperinterval;j++) {
      retvalues[curr_pos]= c_ori_values[skip_c_indices[j]];
      curr_pos++;
    }
    for (int j=0;j<dperinterval;j++) {
      retvalues[curr_pos]= d_ori_values[skip_d_indices[j]];
      curr_pos++;
    }

}
*/

RHSVectors LinearizeKKT::shiftRHSVectors(const RHSVectors& originalRHS,const std::vector<Index>& intervalIDs,const Index& n_intervals, const Index& n_controls,const Index& interval) const
{
  //setup neccessary variables for looping through
  const Index xperinterval = int((originalRHS.static_h_part->Dim()-n_controls)/n_intervals);
  const Index cperinterval = int(originalRHS.static_c_part->Dim()/n_intervals);
  const Index dperinterval = int(originalRHS.static_d_part->Dim()/n_intervals);
  Number* static_h_values = new Number[originalRHS.static_h_part->Dim()-xperinterval];
  Number* u_values = new Number[n_controls];
  Number* static_c_values = new Number[originalRHS.static_c_part->Dim()-cperinterval];
  Number* static_d_values = new Number[originalRHS.static_d_part->Dim()-dperinterval];
  Number* shift_h_values = new Number[xperinterval];
  Number* shift_c_values = new Number[cperinterval];
  Number* shift_d_values = new Number[dperinterval];
  const Number* h_ori_values = originalRHS.static_h_part->Values();
  const Number* c_ori_values = originalRHS.static_c_part->Values();
  const Number* d_ori_values = originalRHS.static_d_part->Values();
  Index* skip_h_indices = new Index[xperinterval];
  Index* skip_c_indices = new Index[cperinterval];
  Index* skip_d_indices = new Index[dperinterval];
  Index* skip_u_indices = new Index[n_controls];
  printf("\nn_controls is: %d\n",n_controls);
  Index curr_pos = 0;
  Index skip_cnt = 0;
  Index u_cnt = 0;

  RHSVectors retval;

  // add objective related data
  printf("\noriginalRHS.static_h_part->Dim()=%d     intervalIDs.size()=%d     interval=%d ,\n",originalRHS.static_h_part->Dim(),intervalIDs.size(),interval);
  for (int j=0;j<originalRHS.static_h_part->Dim();j++) {
    if (intervalIDs[j]==interval) {
      skip_h_indices[skip_cnt]=j;
      skip_cnt++;
    } else if (!intervalIDs[j]) {
      skip_u_indices[u_cnt]=j;
      u_cnt++;
      printf("\nU %d DETECTED!!\n",u_cnt);
    } else {
      static_h_values[curr_pos]=h_ori_values[j];
      curr_pos++;
    }
  }
  curr_pos = 0;
  skip_cnt=0;
  // add equality related data
  for (int j=0;j<originalRHS.static_c_part->Dim();j++) {
    if ((interval-1)*cperinterval<=j && j<interval*cperinterval) {
      skip_c_indices[skip_cnt]=j;
      skip_cnt++;
    }	else {
      static_c_values[curr_pos]=c_ori_values[j];
      curr_pos++;
    }
  }
  curr_pos = 0;
  skip_cnt=0;
  // add inequality related data
  for (int j=0;j<originalRHS.static_d_part->Dim();j++) {
    if ((interval-1)*dperinterval<=j && j<interval*dperinterval) {
      skip_d_indices[skip_cnt]=j;
      skip_cnt++;
    }	else {
      static_d_values[curr_pos]=d_ori_values[j];
      curr_pos++;
    }
  }

  // add the duplicated interval values (objective, equalities, inequalities)
  for (int j=0;j<xperinterval;j++) {
    shift_h_values[j]= h_ori_values[skip_h_indices[j]];
  }
  for (int j=0;j<cperinterval;j++) {
    shift_c_values[j]= c_ori_values[skip_c_indices[j]];
  }
  for (int j=0;j<dperinterval;j++) {
    shift_d_values[j]= d_ori_values[skip_d_indices[j]];
  }

  // add control related data
  for (int j=0;j<n_controls;j++) {
    u_values[j]= h_ori_values[skip_u_indices[j]];
  }

  SmartPtr<DenseVectorSpace> static_h_vec_space = new DenseVectorSpace(originalRHS.static_h_part->Dim()-xperinterval);
  SmartPtr<DenseVectorSpace> static_c_vec_space = new DenseVectorSpace(originalRHS.static_c_part->Dim()-cperinterval);
  SmartPtr<DenseVectorSpace> static_d_vec_space = new DenseVectorSpace(originalRHS.static_d_part->Dim()-dperinterval);
  SmartPtr<DenseVectorSpace> u_vec_space = new DenseVectorSpace(n_controls);
  SmartPtr<DenseVectorSpace> shift_h_vec_space = new DenseVectorSpace(xperinterval);
  SmartPtr<DenseVectorSpace> shift_c_vec_space = new DenseVectorSpace(cperinterval);
  SmartPtr<DenseVectorSpace> shift_d_vec_space = new DenseVectorSpace(dperinterval);

  SmartPtr<DenseVector> static_h_vec = dynamic_cast<const DenseVector*>(static_h_vec_space->MakeNewDenseVector());
  SmartPtr<DenseVector> static_c_vec = dynamic_cast<const DenseVector*>(static_c_vec_space->MakeNewDenseVector());
  SmartPtr<DenseVector> static_d_vec = dynamic_cast<const DenseVector*>(static_d_vec_space->MakeNewDenseVector());
  SmartPtr<DenseVector> u_vec = dynamic_cast<const DenseVector*>(u_vec_space->MakeNewDenseVector());
  SmartPtr<DenseVector> shift_h_vec = dynamic_cast<const DenseVector*>(shift_h_vec_space->MakeNewDenseVector());
  SmartPtr<DenseVector> shift_c_vec = dynamic_cast<const DenseVector*>(shift_c_vec_space->MakeNewDenseVector());
  SmartPtr<DenseVector> shift_d_vec = dynamic_cast<const DenseVector*>(shift_d_vec_space->MakeNewDenseVector());

  static_h_vec->SetValues(static_h_values);
  static_c_vec->SetValues(static_c_values);
  static_d_vec->SetValues(static_d_values);
  u_vec->SetValues(u_values);
  shift_h_vec->SetValues(shift_h_values);
  shift_c_vec->SetValues(shift_c_values);
  shift_d_vec->SetValues(shift_d_values);

  retval.static_h_part = dynamic_cast<const DenseVector*>(GetRawPtr(shift_h_vec));
  retval.static_c_part = dynamic_cast<const DenseVector*>(GetRawPtr(shift_c_vec));
  retval.static_d_part = dynamic_cast<const DenseVector*>(GetRawPtr(shift_d_vec));
  retval.u_part = dynamic_cast<const DenseVector*>(GetRawPtr(u_vec));
  retval.shift_h_part = dynamic_cast<const DenseVector*>(GetRawPtr(shift_h_vec));
  retval.shift_c_part = dynamic_cast<const DenseVector*>(GetRawPtr(shift_c_vec));
  retval.shift_d_part = dynamic_cast<const DenseVector*>(GetRawPtr(shift_d_vec));


  return retval;
}

RHSVectors LinearizeKKT::getShiftedRHS(const Index& column,SmartPtr<IpoptApplication> app,const Index & interval) const
{
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  // get raw material for multicolumn matrix rhs
  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<const Vector> y_c = app->IpoptDataObject()->curr()->y_c();
  SmartPtr<const Vector> y_d = app->IpoptDataObject()->curr()->y_d();
  SmartPtr<const Matrix> rhs_c = orig_nlp->jac_c_p(*x);
  SmartPtr<const Matrix> rhs_d = orig_nlp->jac_d_p(*x);
  SmartPtr<const Matrix> rhs_h = orig_nlp->h_p(*x, 1.0, *y_c, *y_d);

  SmartPtr<const VectorSpace> x_space = x->OwnerSpace();
  SmartPtr<const MatrixSpace> rhs_h_space = rhs_h->OwnerSpace();
  SmartPtr<MultiVectorMatrix> h_MVM =  getMVMFromMatrix(rhs_h);
  SmartPtr<MultiVectorMatrix> c_MVM =  getMVMFromMatrix(rhs_c);
  SmartPtr<MultiVectorMatrix> d_MVM =  getMVMFromMatrix(rhs_d);

  const IntervalInfoSet intervals = IntervalInfoSet(dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p())));

  const std::vector<Index> intervalIDs = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_space))->GetIntegerMetaData("intervalID");

  //  printf("\n intervalIDs.size()=%d\n",intervalIDs.size());
  //  for (int i=0;i<intervalIDs.size();i++)
  //    printf("intervalIDs[%d]=%d\n",i,intervalIDs[i]);
  const Index n_intervals = intervals.getIntervalCount();
  const Index xperinterval = getEntryCountPerInterval(x);
  const Index n_u = x->Dim()-n_intervals*xperinterval;
  assert(!column<orig_nlp->p()->Dim()||column<0);
  RHSVectors original;
  original.static_h_part = dynamic_cast<const DenseVector*>(GetRawPtr(h_MVM->GetVector(column)));
  original.u_part = NULL;
  original.static_c_part = dynamic_cast<const DenseVector*>(GetRawPtr(c_MVM->GetVector(column)));
  original.static_d_part = dynamic_cast<const DenseVector*>(GetRawPtr(d_MVM->GetVector(column)));
  original.shift_h_part =NULL;
  original.shift_c_part =NULL;
  original.shift_d_part =NULL;
  return this->shiftRHSVectors(original,intervalIDs,n_intervals,n_u,interval);

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  // initialize rhs_matrix
  //  const IntervalInfoSet intervals = IntervalInfoSet(dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p())));
  const Index n_p = orig_nlp->p()->Dim();
  //  const Index n_intervals = intervals.getIntervalCount();
  // const Index xperinterval = getEntryCountPerInterval(x);
  // const Index n_u = x->Dim()-n_intervals*xperinterval;
  const Index cperinterval = int(rhs_c->NRows()/n_intervals);
  const Index dperinterval = int(rhs_d->NRows()/n_intervals);
  const Index n_rows = int(x->Dim()+xperinterval+(cperinterval+dperinterval)*(n_intervals+1));
  SmartPtr<const DenseVectorSpace> retvec_space = new DenseVectorSpace(n_rows);
  SmartPtr<MultiVectorMatrixSpace> retval_space = new MultiVectorMatrixSpace(n_p,*retvec_space);
  SmartPtr<MultiVectorMatrix> retval = retval_space->MakeNewMultiVectorMatrix();

  // setup quantities neccessary for rhs setup


  Index curr_pos = 0;
  SmartPtr<const DenseVector> h_ori_vec;
  SmartPtr<const DenseVector> c_ori_vec;
  SmartPtr<const DenseVector> d_ori_vec;
  Number* retvalues = new Number[n_rows];
  Index* skip_h_indices = new Index[xperinterval];
  Index* skip_c_indices = new Index[cperinterval];
  Index* skip_d_indices = new Index[dperinterval];
  printf("\nn_u is: %d\n",n_u);
  Index* skip_u_indices = new Index[n_u];
  Index skip_cnt = 0;
  Index u_cnt = 0;

  // fill rhs_matrix with values
  for (int i=0;i<n_p;i++) {
    SmartPtr<DenseVector> retvec = dynamic_cast<const DenseVector*>(retvec_space->MakeNewDenseVector());
    h_ori_vec = dynamic_cast<const DenseVector*>(GetRawPtr(h_MVM->GetVector(i)));
    const Number* h_ori_values = h_ori_vec->Values();
    c_ori_vec = dynamic_cast<const DenseVector*>(GetRawPtr(c_MVM->GetVector(i)));
    const Number* c_ori_values = c_ori_vec->Values();
    d_ori_vec = dynamic_cast<const DenseVector*>(GetRawPtr(d_MVM->GetVector(i)));
    const Number* d_ori_values = d_ori_vec->Values();

    //    h_ori_vec->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "h_ori_vec");
    c_ori_vec->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "c_ori_vec");
    d_ori_vec->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "d_ori_vec");
    // add objective related data
    for (int j=0;j<x->Dim();j++) {
      if (intervalIDs[j]==interval) {
	skip_h_indices[skip_cnt]=j;
	skip_cnt++;
      } else if (!intervalIDs[j]) {
	skip_u_indices[u_cnt]=j;
	u_cnt++;
	printf("\nU %d DETECTED!!\n",u_cnt);
      } else {
	retvalues[curr_pos]=h_ori_values[j];
	curr_pos++;
      }
    }
    skip_cnt=0;
    // add equality related data
    for (int j=0;j<rhs_c->NRows();j++) {
      if ((interval-1)*cperinterval<=j && j<interval*cperinterval) {
	skip_c_indices[skip_cnt]=j;
	skip_cnt++;
      }	else {
	retvalues[curr_pos]=c_ori_values[j];
	curr_pos++;
      }
    }
    skip_cnt=0;
    // add inequality related data
    for (int j=0;j<rhs_d->NRows();j++) {
      if ((interval-1)*dperinterval<=j && j<interval*dperinterval) {
	skip_d_indices[skip_cnt]=j;
	skip_cnt++;
      }	else {
	retvalues[curr_pos]=d_ori_values[j];
	curr_pos++;
      }
    }

    // add the duplicated interval values (objective, equalities, inequalities)
    for (int j=0;j<xperinterval;j++) {
      retvalues[curr_pos]= h_ori_values[skip_h_indices[j]];
      curr_pos++;
    }
    for (int j=0;j<cperinterval;j++) {
      retvalues[curr_pos]= c_ori_values[skip_c_indices[j]];
      curr_pos++;
    }
    for (int j=0;j<dperinterval;j++) {
      retvalues[curr_pos]= d_ori_values[skip_d_indices[j]];
      curr_pos++;
    }

    // add control related data
    for (int j=0;j<n_u;j++) {
      printf("\nadding control related data: curr_pos = %d     skip_u_indices[%d]=%d    h_ori_values.size()=%d ",curr_pos,j,skip_u_indices[j],h_ori_vec->Dim());
      retvalues[curr_pos]= h_ori_values[skip_u_indices[j]];
      curr_pos++;
    }

    // add the duplicated interval values (objective, equalities, inequalities)
    for (int j=0;j<xperinterval;j++) {
      retvalues[curr_pos]= h_ori_values[skip_h_indices[j]];
      curr_pos++;
    }
    for (int j=0;j<cperinterval;j++) {
      retvalues[curr_pos]= c_ori_values[skip_c_indices[j]];
      curr_pos++;
    }
    for (int j=0;j<dperinterval;j++) {
      retvalues[curr_pos]= d_ori_values[skip_d_indices[j]];
      curr_pos++;
    }

    retval->SetVector(i,*retvec);
    curr_pos = 0;
    skip_cnt = 0;
    u_cnt = 0;
 }


}

SmartPtr<Matrix> LinearizeKKT::getShiftedSensitivities(SmartPtr<IpoptApplication> app,const Index & interval) const
{
  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  IntervalInfoSet intervals = IntervalInfoSet(p);
  SmartPtr<OptionsList> options = app->Options();
  Index np = orig_nlp->p()->Dim();
  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<MultiVectorMatrixSpace> rhs_space = new MultiVectorMatrixSpace(np, *x->OwnerSpace());
  SmartPtr<MultiVectorMatrix> mv = dynamic_cast<MultiVectorMatrix*>(rhs_space->MakeNew());

  return sens_matrix;
}

ParameterShift* assignShiftMethod(SmartPtr<OptionsList> options)
{
  // no other shift methods implemented yet!
  return new LinearizeKKT();
}

////////////////////////////end of intervallization pseudo implementation///////////////////////////


int main(int argc, char**args)
{

  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Check if executable run only to print out options documentation
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
  // Add the suffix for interval-organization
  suffix_handler->AddAvailableSuffix("intervalID", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);

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
  //opt_jac_c_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_c_p");
  //opt_jac_d_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_d_p");
  //opt_h_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_h_p");

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
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  IntervalInfoSet intervals = IntervalInfoSet(p);
  SmartPtr<OptionsList> options = app->Options();
  BranchingCriterion* branchmode = assignBranchingMethod(options);
  std::vector<SplitChoice> splitchoices = branchmode->branchSensitivityMatrix(app);

  ControlSelector* pickfirst = assignControlMethod(options);
  SplitDecision resulting_split = pickfirst->decideSplitControl(splitchoices);
  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  ParameterShift* shifter = assignShiftMethod(options);
  RHSVectors shifted_rhs = (shifter)->getShiftedRHS(0,app,1);

  printf("\nproposed for split: intervalID: %d parameter: %d\n",resulting_split.intervalID,resulting_split.parameterID);

  //write gathered information into .dat file to access with python
  std::string fname = "branch_intervals.dat";
  std::ofstream branch_intervals;
  branch_intervals.open(fname.c_str());
  char buffer[63];
  branch_intervals << "#.dat file automatically generated by AMPL intervallization routine\n#Ben Waldecker Sep 2012\n";
  sprintf(buffer,"\nintervalID: %d parameter: %d\n",resulting_split.intervalID,resulting_split.parameterID);
  branch_intervals << buffer;
  branch_intervals << "\n\n#end of file";
  branch_intervals.close();

  return 1;
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

