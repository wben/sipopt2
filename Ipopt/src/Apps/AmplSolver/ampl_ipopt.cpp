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
  IntervalInfo(const Number value, const Index parameterID, const Index intervalID, const Index vector_index, const bool is_upper);
  ~IntervalInfo();
  void setParameters(const std::vector<std::string> pnames, const std::vector<Number> pvalues);
  void addParameter(const std::vector<std::string> pnames, const std::vector<Number> pvalues);
  void getIndex(Index &index);
  void getValue(Number &value);
  void setValue(Number value);
  void getIntervalID (Index &nint);
  void getParameterID (Index &paraID);
  bool isUpper();
  void setIntervals(const Index nint);
  void printSet();

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

  IntervalInfoSet(std::vector<IntervalInfo> intinfovec);
  IntervalInfoSet(SmartPtr<const DenseVector> parameters);
  ~IntervalInfoSet();
  void setIntInfoSet(std::vector<IntervalInfo> intinfovec);
  void getIntInfoSet(std::vector<IntervalInfo> &intinfovec);
  void getIndexVec(std::vector<Index> &indexvec);
  void getIndex(Index intindex, Index &index);
  void getValueVec(std::vector<Number> &valuevec);
  void getValue(Index intindex, Number &value);
  void setValueVec(std::vector<Number> valuevec);
  void setValue(Index intindex, Number value);
  void getIntervalIDVec (std::vector<Index> &nintvec);
  void getIntervalID(Index intindex, Index &intervalID);
  void getIntervalIDs(std::vector<Index>&intervalIDs);
  void getParameterIDVec (std::vector<Index> &parameterIDvec);
  void getParameterID(Index paraindex, Index &parameterID);
  void getParameterIDs(std::vector<Index> &parameterIDs);
  void isUpperVec(std::vector<bool> &is_uppervec);
  bool isUpper(Index isupperindex);
  void getOtherBndIdx(Index boundindex,Index &otherbndidx);
  void getParameterCount(Index &paracount);
  void getIntervalCount(Index &intervalcount);
  void printSet();
  Index Size();

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
  virtual std::vector<Number> scaleIntervalWidths(IntervalInfoSet intervals,Index intervalID) = 0;
};
IntervalWidthScaling* assignScalingMethod(SmartPtr<OptionsList> options);

class NoScaling : public IntervalWidthScaling
{
public:
  virtual std::vector<Number> scaleIntervalWidths(IntervalInfoSet intervals,Index intervalID);
};

class TotIntWidthScaling : public IntervalWidthScaling
{
public:
  virtual std::vector<Number> scaleIntervalWidths(IntervalInfoSet intervals,Index intervalID);
};

class IntWidthScaling : public IntervalWidthScaling
{
public:
  virtual std::vector<Number> scaleIntervalWidths(IntervalInfoSet intervals,Index intervalID);
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

class BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app) = 0;
  virtual SplitChoice branchInterval(std::vector<SplitChoice> splitchoices) = 0;
};

BranchingCriterion* assignBranchingMethod(SmartPtr<OptionsList> options);

class RandomBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app);
  virtual SplitChoice branchInterval(std::vector<SplitChoice> splitchoices);
};

class LargerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app);
  virtual SplitChoice branchInterval(std::vector<SplitChoice> splitchoices);
};

class SmallerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app);
  virtual SplitChoice branchInterval(std::vector<SplitChoice> splitchoices);
};

class AbsoluteLargerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app);
  virtual SplitChoice branchInterval(std::vector<SplitChoice> splitchoices);
};

class Intervaluation
{
public:
  virtual SplitChoice intervaluateInterval(Index controlrow, Index intervalID,SmartPtr<IpoptApplication> app) = 0;
};
Intervaluation* assignIntervaluationMethod(SmartPtr<OptionsList> options);

class OneBoundIntervaluation : public Intervaluation
{
public:
  SplitChoice intervaluateInterval(Index controlrow, Index intervalID,SmartPtr<IpoptApplication> app);
};

class BothBoundIntervaluation : public Intervaluation
{
public:
  SplitChoice intervaluateInterval(Index controlrow, Index intervalID,SmartPtr<IpoptApplication> app);
};

class ControlSelector
{
public:
  virtual SplitDecision decideSplitControl(std::vector<SplitChoice> choices) = 0;
};

class SelectFirstControl : public ControlSelector
{
public:
  SplitDecision decideSplitControl(std::vector<SplitChoice> choices);
};

ControlSelector* assignControlMethod(SmartPtr<OptionsList> options);

/////////////////////////////end of intervallization pseudo header///////////////////////////////////

SmartPtr<Matrix> getSensitivityMatrix(SmartPtr<IpoptApplication> app);
SmartPtr<Vector> getDirectionalDerivative(SmartPtr<IpoptApplication> app,
					  SmartPtr<Matrix> sens_matrix);
bool doIntervalization(SmartPtr<IpoptApplication> app);
std::vector<Number> getIntervalWidths(IntervalInfoSet intervals, bool do_scaling=false);
Number Abs(Number value);

///////////////////////////start of intervallization pseudo implementation//////////////////////////
IntervalInfo::IntervalInfo() {}

IntervalInfo::IntervalInfo(const Number value, const Index parameterID, const Index intervalID, const Index vector_index, const bool is_upper)
{

  value_ = value;
  parameterID_ = parameterID;
  intervalID_ = intervalID;
  index_ = vector_index;
  is_upper_ = is_upper;

}

IntervalInfo:: ~IntervalInfo() {}

void IntervalInfo::setParameters(const std::vector<std::string> pnames, const std::vector<Number> pvalues)
{  }

void IntervalInfo::addParameter(const std::vector<std::string> pnames, const std::vector<Number> pvalues)
{  }

void IntervalInfo::getIndex(Index &index)
{
  index = index_;
}

void IntervalInfo::getValue(Number &value)
{
  value = value_;
}

void IntervalInfo::setValue(Number value)
{
  value_=value;
}

void IntervalInfo::getIntervalID(Index &nint)
{
  nint = intervalID_;
}

void IntervalInfo::getParameterID(Index &paraID)
{
  paraID = parameterID_;
}

bool IntervalInfo::isUpper()
{
  return is_upper_;
}

void IntervalInfo::setIntervals(const Index nint)
{  }

void IntervalInfo::printSet()
{
  printf("\n value: %f parameterID: %d, intervalID: %d, index: %d, is_upper: %d \n", value_, parameterID_, intervalID_, index_, is_upper_);
}


IntervalInfoSet::IntervalInfoSet() { }

IntervalInfoSet::IntervalInfoSet(std::vector<IntervalInfo> intinfovec)
{
  valuevec_.clear();
  intinfovec_.clear();
  indexvec_.clear();
  parameterIDvec_.clear();
  intervalIDvec_.clear();
  is_uppervec_.clear();
  Index tmp_index;
  Index tmp_paraID;
  Index tmp_intID;
  int i=0;

  // sort algorithm to make sure entry indexing in IntervalInfoSet matches given indexing
  while (intinfovec_.size()<intinfovec.size()) {
    intinfovec[i].getIndex(tmp_index);
    if (tmp_index==indexvec_.size()) {
      intinfovec_.push_back(intinfovec[i]);
      indexvec_.push_back(tmp_index);
      intinfovec[i].getParameterID(tmp_paraID);
      parameterIDvec_.push_back(tmp_paraID);
      intinfovec[i].getIntervalID(tmp_intID);
      intervalIDvec_.push_back(tmp_intID);
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

void IntervalInfoSet::setIntInfoSet(std::vector<IntervalInfo> intinfovec)
{
  intinfovec_=intinfovec;
}

void IntervalInfoSet::getIntInfoSet(std::vector<IntervalInfo> &intinfovec)
{
  intinfovec=intinfovec_;
}

void IntervalInfoSet::getIndexVec(std::vector<Index> &indexvec)
{
  indexvec=indexvec_;
}

void IntervalInfoSet::getIndex(Index intindex, Index &index)
{
  if (intindex<indexvec_.size())
    index = indexvec_[intindex];
  else {
    index = -1;
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::getIndex() call with out of range index!\n");
  }
}

void IntervalInfoSet::getValueVec(std::vector<Number> &valuevec)
{
  valuevec = valuevec_;
}

void IntervalInfoSet::getValue(Index intindex,Number &value)
{
  intinfovec_[intindex].getValue(value);
}

void IntervalInfoSet::setValueVec(std::vector<Number> valuevec)
{
  valuevec_ = valuevec;
}

void IntervalInfoSet::setValue(Index intindex,Number value)
{
  intinfovec_[intindex].setValue(value);
}

void IntervalInfoSet::getIntervalIDVec(std::vector<Index> &intervalIDvec)
{
  intervalIDvec = intervalIDvec_;
}

void IntervalInfoSet::getIntervalID(Index intindex, Index &intervalID)
{
  if (intindex<intervalIDvec_.size())
    intervalID = intervalIDvec_[intindex];
  else {
    intervalID = -1;
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::getIntervalID() call with out of range index!\n");
  }
}

void IntervalInfoSet::getIntervalIDs(std::vector<Index>&intervalIDs)
{
  std::vector<Index> tmp_IDs;
  Index nints;
  this->getIntervalCount(nints);
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
  intervalIDs = tmp_IDs;
}

void IntervalInfoSet::getParameterIDVec(std::vector<Index> &parameterIDvec)
{
  parameterIDvec = parameterIDvec_;
}

void IntervalInfoSet::getParameterID(Index paraindex, Index &parameterID)
{
  if (paraindex<parameterIDvec_.size())
    parameterID = parameterIDvec_[paraindex];
  else {
    parameterID = -1;
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::getParameterID() call with out of range index!\n");
  }
}

void IntervalInfoSet::getParameterIDs(std::vector<Index>&parameterIDs)
{
  std::vector<Index> tmp_IDs;
  Index nints;
  this->getParameterCount(nints);
  Index* tmp_new_ID = new Index;
  Index* tmp_ID = new Index;
  *tmp_new_ID = parameterIDvec_[0];
  // assign first value
  for (int i=0;i<parameterIDvec_.size();i++)
    if (parameterIDvec_[i]<*tmp_new_ID)
      *tmp_new_ID=parameterIDvec_[i];
  tmp_IDs.push_back(*tmp_new_ID);

  //assign follow up values
  while (tmp_IDs.size()<nints) {
    *tmp_ID = *tmp_new_ID;
    for (int i=0;i<parameterIDvec_.size();i++) {
      if (parameterIDvec_[i]>*tmp_ID && (*tmp_ID==*tmp_new_ID))
 	*tmp_new_ID=parameterIDvec_[i];
      if (parameterIDvec_[i]>*tmp_ID && parameterIDvec_[i]<*tmp_new_ID)
	*tmp_new_ID = parameterIDvec_[i];
    }
    tmp_IDs.push_back(*tmp_new_ID);
  }
  parameterIDs = tmp_IDs;
}

void IntervalInfoSet::isUpperVec(std::vector<bool> &is_uppervec)
{
  is_uppervec=is_uppervec_;
}

bool IntervalInfoSet::isUpper(Index isupperindex)
{
  if (isupperindex<is_uppervec_.size())
    return is_uppervec_[isupperindex];
  else {
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::isUpper() call with out of range index!\n");
    printf("\nsize of is_uppervec_: %d \n", is_uppervec_.size());
    return 0;
  }
}
void IntervalInfoSet::getOtherBndIdx(Index boundindex,Index &otherbndidx)
{
  for (int i=0;i<intinfovec_.size();i++){
    if (parameterIDvec_[int(boundindex)]==parameterIDvec_[i] && intervalIDvec_[int(boundindex)]==intervalIDvec_[i] && is_uppervec_[int(boundindex)]!=is_uppervec_[i]){
      otherbndidx = indexvec_[i];
      i=intinfovec_.size();
    }
  }
}

void IntervalInfoSet::getParameterCount(Index &paracount)
{
  Index tmp_count =0;
  for (int i=0;i<parameterIDvec_.size();i++) {
    if (i==0)
      tmp_count = parameterIDvec_[i];
    if (parameterIDvec_[i]>tmp_count)
      tmp_count = parameterIDvec_[i];
  }
  paracount = tmp_count;
}

void IntervalInfoSet::getIntervalCount(Index &intervalcount)
{
  Index tmp_count= 0;
  for (int i=0;i<intervalIDvec_.size();i++) {
    if (i==0)
      tmp_count = intervalIDvec_[i];
    if (intervalIDvec_[i]>tmp_count)
      tmp_count = intervalIDvec_[i];
  }
  intervalcount = tmp_count;
}

void IntervalInfoSet::printSet()
{
  for (int i=0; i<intinfovec_.size();i++){
    printf("\n\nIntervalInfoSet Eintrag %d:\n", i);
    intinfovec_[i].printSet();
    printf("\n");
  }
}

Index IntervalInfoSet::Size()
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


std::vector<Number> NoScaling::scaleIntervalWidths(IntervalInfoSet intervals,Index intervalID)
{
  // no scaling: just return a vector with apropriate length containing 1s
  Index npars;
  intervals.getParameterCount(npars);
  std::vector<Number> retval(npars);
  for (int i=0;i<npars;i++)
    retval[i]=1;
  return retval;
}

std::vector<Number> TotIntWidthScaling::scaleIntervalWidths(IntervalInfoSet intervals,Index intervalID)
{
  Index npars;
  intervals.getParameterCount(npars);
  std::vector<Number> retval(npars);
  for (int i=0;i<npars;i++)
    retval[i]=0;
  std::vector<Number> par_values;
  intervals.getValueVec(par_values);
  std::vector<Index> intervalflags;
  intervals.getIntervalIDVec(intervalflags);
  std::vector<Index> parameterflags;
  intervals.getParameterIDVec(parameterflags);
  std::vector<Index> parameterIDs;
  intervals.getParameterIDs(parameterIDs);
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

std::vector<Number> IntWidthScaling::scaleIntervalWidths(IntervalInfoSet intervals,Index intervalID)
{
  Index npars;
  intervals.getParameterCount(npars);
  std::vector<Number> retval(npars);
  for (int i=0;i<npars;i++)
    retval[i]=0;
  std::vector<Number> par_values;
  intervals.getValueVec(par_values);
  std::vector<Index> intervalflags;
  intervals.getIntervalIDVec(intervalflags);
  std::vector<Index> parameterflags;
  intervals.getParameterIDVec(parameterflags);
  std::vector<Index> parameterIDs;
  intervals.getParameterIDs(parameterIDs);

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

std::vector<SplitChoice> RandomBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app)
{
  // RandomBranching not implemented yet.
}

SplitChoice RandomBranching::branchInterval(std::vector<SplitChoice> splitchoices)
{
  // RandomBranching not implemented yet.
}

std::vector<SplitChoice>  LargerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app)
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

  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs;
  intervals.getIntervalIDs(intervalIDs);
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater = assignIntervaluationMethod(options);
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

SplitChoice LargerBranching::branchInterval(std::vector<SplitChoice> splitchoices)
{
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

std::vector<SplitChoice> SmallerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app)
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

  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs;
  intervals.getIntervalIDs(intervalIDs);
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater = assignIntervaluationMethod(options);
  for (Index i=0;i<ctrl_rows.size();i++) {
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

SplitChoice SmallerBranching::branchInterval(std::vector<SplitChoice> splitchoices)
{
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

std::vector<SplitChoice> AbsoluteLargerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app)
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

  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs;
  intervals.getIntervalIDs(intervalIDs);
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater = assignIntervaluationMethod(options);
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

SplitChoice AbsoluteLargerBranching::branchInterval(std::vector<SplitChoice> splitchoices)
{
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

SplitChoice OneBoundIntervaluation::intervaluateInterval(Index controlrow, Index intervalID,SmartPtr<IpoptApplication> app)
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
  std::vector<Number> p_values;
  intervals.getValueVec(p_values);
  std::vector<Index> ID_vec;
  intervals.getIntervalIDVec(ID_vec);
  std::vector<Index> para_vec;
  intervals.getParameterIDVec(para_vec);
  Index npars;
  intervals.getParameterCount(npars);
  std::vector<SplitChoice> int_splitchoices(2*npars);
  std::vector<Index> parameterIDs;
  intervals.getParameterIDs(parameterIDs);
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
	intervals.getOtherBndIdx(j,tmp_idx);
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

SplitChoice BothBoundIntervaluation::intervaluateInterval(Index controlrow, Index intervalID,SmartPtr<IpoptApplication> app)
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
  std::vector<Number> p_values;
  intervals.getValueVec(p_values);
  std::vector<Index> ID_vec;
  intervals.getIntervalIDVec(ID_vec);
  std::vector<Index> para_vec;
  intervals.getParameterIDVec(para_vec);
  Index npars;
  intervals.getParameterCount(npars);
  std::vector<SplitChoice> int_splitchoices(npars);
  std::vector<Index> parameterIDs;
  intervals.getParameterIDs(parameterIDs);
  Number tmp_value;
  Index tmp_idx;
  //cycle through the different parameter values for the given interval and determine branching objective values
  for (int i=0; i<npars;i++) {
    int_splitchoices[i].parameterID = parameterIDs[i];
    int_splitchoices[i].intervalID = intervalID;
    for (int j=0;j<intervals.Size();j++) {
      if (ID_vec[j]==intervalID && para_vec[j]==parameterIDs[i]) {
      // parameter entry for interval in question found!
	intervals.getOtherBndIdx(j,tmp_idx);
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

SplitDecision SelectFirstControl::decideSplitControl(std::vector<SplitChoice> choices)
{
  //simply returns first entry of SplitDecisions
  SplitDecision retval;
  retval.intervalID = choices[0].intervalID;
  retval.parameterID = choices[0].parameterID;
  return retval;
}

ControlSelector* assignControlMethod(SmartPtr<OptionsList> options)
{
  return new SelectFirstControl();
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
  //  opt_jac_c_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_c_p");
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

