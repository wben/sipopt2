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
  virtual void scaleIntervalWidths(IntervalInfoSet intervals,std::vector<Number>& intervalwidths) = 0;
};
IntervalWidthScaling* assignScalingMethod(SmartPtr<OptionsList> options);

class NoScaling : public IntervalWidthScaling
{
public:
  virtual void scaleIntervalWidths(IntervalInfoSet intervals, std::vector<Number>& intervalwidths);
};

class TotIntWidthScaling : public IntervalWidthScaling
{
public:
  virtual void scaleIntervalWidths(IntervalInfoSet intervals, std::vector<Number>& intervalwidths);
};

class IntWidthScaling : public IntervalWidthScaling
{
public:
  virtual void scaleIntervalWidths(IntervalInfoSet intervals, std::vector<Number>& intervalwidths);
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
  virtual SplitChoice branchInterval(Index intervalID,SmartPtr<IpoptApplication> app) = 0;
};

BranchingCriterion* assignBranchingMethod(SmartPtr<OptionsList> options);

class RandomBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app);
  virtual SplitChoice branchInterval(Index intervalID,SmartPtr<IpoptApplication> app);
};

class LargerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app);
  virtual SplitChoice branchInterval(Index intervalID,SmartPtr<IpoptApplication> app);
};

class SmallerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app);
  virtual SplitChoice branchInterval(Index intervalID,SmartPtr<IpoptApplication> app);
};

class AbsoluteLargerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app);
  virtual SplitChoice branchInterval(Index intervalID,SmartPtr<IpoptApplication> app);
};

class Intervaluation
{
public:
  virtual SplitChoice intervaluateInterval(Index intervalID,SmartPtr<IpoptApplication> app) = 0;
};
Intervaluation* assignIntervaluationMethod(SmartPtr<OptionsList> options);

class OneBoundIntervaluation : public Intervaluation
{
public:
  SplitChoice intervaluateInterval(Index intervalID,SmartPtr<IpoptApplication> app);
};

class BothBoundIntervaluation : public Intervaluation
{
public:
  SplitChoice intervaluateInterval(Index intervalID,SmartPtr<IpoptApplication> app);
};

class ControlSelector
{
public:
  virtual SplitDecision decideSplitControl(std::vector<SplitChoice> choices) = 0;
};

class SelectFirstControl : public ControlSelector
{
public:
  virtual SplitDecision decideSplitControl(std::vector<SplitChoice> choices);
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
  Index* tmp_new_ID = new Index;
  Index* tmp_ID = new Index;
  *tmp_ID = intervalIDvec_[0];
  *tmp_new_ID = intervalIDvec_[0];
  // assign first value
  for (int i=0;i<intervalIDvec_.size();i++) {
    if (intervalIDvec_[i]<*tmp_new_ID)
      *tmp_new_ID=intervalIDvec_[i];
  }
  tmp_IDs.push_back(*tmp_new_ID);
  *tmp_ID = tmp_IDs[0];
  *tmp_new_ID = *tmp_ID +1;
  while (*tmp_ID!=*tmp_new_ID) {
    for (int i=0;i<intervalIDvec_.size();i++) {
      if (intervalIDvec_[i]>*tmp_ID && intervalIDvec_[i]<*tmp_new_ID)
	*tmp_new_ID=intervalIDvec_[i];
    }
    tmp_IDs.push_back(*tmp_new_ID);
    *tmp_ID = *tmp_new_ID;
  }

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
    printf("\nscalingmode is set to %s\n",scalingmode.c_str());
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

void NoScaling::scaleIntervalWidths(IntervalInfoSet intervals, std::vector<Number>& intervalwidths)
{
  std::vector<Number> tmp_widths(intervals.Size());
  for (int i=0;i<tmp_widths.size();i++)
    tmp_widths[i]=1;
  intervalwidths=tmp_widths;
}

void TotIntWidthScaling::scaleIntervalWidths(IntervalInfoSet intervals, std::vector<Number>& intervalwidths)
{

  std::vector<Number> par_values;
  intervals.getValueVec(par_values);

  std::vector<Number> tmp_widths(intervals.Size());

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
  std::vector<Index> lower_idx(intervals.Size());
  std::vector<Index> upper_idx(intervals.Size());
  int tmp_j =0;
  for (Index i=0;i<intervals.Size();i++) {
    if (intervals.isUpper(i)) {
      intervals.getIndex(i,upper_idx[tmp_j]);
      intervals.getOtherBndIdx(upper_idx[tmp_j],lower_idx[tmp_j]);
      tmp_j++;
    } else {
      intervals.getIndex(i,lower_idx[tmp_j]);
      intervals.getOtherBndIdx(lower_idx[tmp_j],upper_idx[tmp_j]);
      tmp_j++;
    }
  }

  // determine intervalwidths for each parametervalue - the index of the intervalwidth stored here matches the one of the parameterdata in IntervalInfoSet intervals
  Number tmp_width;
  for (Index i=0;i<tmp_widths.size();i++) {
    intervals.getParameterID(upper_idx[i],*tmp_par);
    tmp_width=(par_values[upper_idx.at(i)]-par_values[lower_idx.at(i)])/total_int_widths[*tmp_par-1];
    tmp_widths[upper_idx.at(i)] = tmp_width;
    tmp_widths[lower_idx.at(i)] = tmp_width;
    // printf("\n skalierte Intervalbreite %d für Parameter %d beträgt: %f\n",i,*tmp_par,intervalwidths[i]);
  }
  for (int i=0;i<tmp_widths.size();i++) {
    printf("\ntmp_widths.size ist: %d. upper_idx.size() ist: %d. lower_idx.size() ist: %d.\n",tmp_widths.size(),upper_idx.size(),lower_idx.size());
    printf("\ntmp_widths.[%d] ist: %f\n", i,tmp_widths[i]);
    printf("\nupper_idx.[%d] ist: %d\n", i,upper_idx[i]);
    printf("\nlower_idx.[%d] ist: %d\n", i,lower_idx[i]);
  }
  intervalwidths = tmp_widths;
}

void IntWidthScaling::scaleIntervalWidths(IntervalInfoSet intervals, std::vector<Number>& intervalwidths)
{

  std::vector<Number> par_values;
  intervals.getValueVec(par_values);

  std::vector<Number> tmp_widths(intervals.Size());

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
  //cycle through intervalset to assign upper and lower value vector indexes to each other
  std::vector<Index> lower_idx(intervals.Size());
  std::vector<Index> upper_idx(intervals.Size());
  int tmp_j =0;
  for (Index i=0;i<intervals.Size();i++) {
    if (intervals.isUpper(i)) {
      intervals.getIndex(i,upper_idx[tmp_j]);
      intervals.getOtherBndIdx(upper_idx[tmp_j],lower_idx[tmp_j]);
      tmp_j++;
    } else {
      intervals.getIndex(i,lower_idx[tmp_j]);
      intervals.getOtherBndIdx(lower_idx[tmp_j],upper_idx[tmp_j]);
      tmp_j++;
    }
  }

  // determine intervalwidths for each parametervalue - the index of the intervalwidth stored here matches the one of the parameterdata in IntervalInfoSet intervals
  Number tmp_width;
  for (Index i=0;i<tmp_widths.size();i++) {
    intervals.getParameterID(upper_idx[i],*tmp_par);
    tmp_width=(par_values[upper_idx.at(i)]-par_values[lower_idx.at(i)]);
    tmp_widths[upper_idx.at(i)] = tmp_width;
    tmp_widths[lower_idx.at(i)] = tmp_width;
    // printf("\n skalierte Intervalbreite %d für Parameter %d beträgt: %f\n",i,*tmp_par,intervalwidths[i]);
  }
  for (int i=0;i<tmp_widths.size();i++) {
    printf("\ntmp_widths.size ist: %d. upper_idx.size() ist: %d. lower_idx.size() ist: %d.\n",tmp_widths.size(),upper_idx.size(),lower_idx.size());
    printf("\ntmp_widths.[%d] ist: %f\n", i,tmp_widths[i]);
    printf("\nupper_idx.[%d] ist: %d\n", i,upper_idx[i]);
    printf("\nlower_idx.[%d] ist: %d\n", i,lower_idx[i]);
  }
  intervalwidths = tmp_widths;
}

BranchingCriterion* assignBranchingMethod(SmartPtr<OptionsList> options)
{
  std::string branchmode;
  if (options->GetStringValue("branchmode",branchmode ,"")){
    printf("\branchmode is set to %s\n",branchmode.c_str());
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

}

SplitChoice RandomBranching::branchInterval(Index intervalID, SmartPtr<IpoptApplication> app)
{

}

std::vector<SplitChoice>  LargerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app)
{

}

SplitChoice LargerBranching::branchInterval(Index intervalID, SmartPtr<IpoptApplication> app)
{

}

std::vector<SplitChoice> SmallerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app)
{
  SmartPtr<OptionsList> options = app->Options();
  Intervaluation* intervaluate = assignIntervaluationMethod(options);
  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));

  const Index nrows = mv_sens->NRows();
  const Index ncols = mv_sens->NCols();
  // vector to store the indexes of sensitivity matrix rows containing control related data
  std::vector<Index> ctrl_rows;


  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  SmartPtr<const DenseVectorSpace> p_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(p->OwnerSpace()));

  /*  // get parameter names
  const std::vector<std::string> parnames = p_space->GetStringMetaData("idx_names");
  const Index i_p = p_space->Dim();
  std::vector<std::string> par_names_tmp;
  for (int i=0;i<i_p;i++)
    par_names_tmp.push_back(parnames[i].c_str());
  const std::vector<std::string> par_names = par_names_tmp;
  */

  IntervalInfoSet intervals = IntervalInfoSet(p);
  // get parameter values
  std::vector<Number> par_values;
  intervals.getValueVec(par_values)
  std::vector<Number> intervalwidths;
  IntervalWidthScaling* scalingmode = assignScalingMethod(options);
  scalingmode->scaleIntervalWidths(intervals,intervalwidths);
  for (int i=0; i<intervalwidths.size();i++)
    printf("\n skalierte Intervalbreite %d beträgt: %f\n",i,intervalwidths[i]);

// cycle through var space interval flags to identify and save control indexes
  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");
  for (int i=0;i<nrows;i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
    //  // printf("\ndx/dp MetaData: intervalID an der Stelle %d hat den Wert %d.\n", i, var_int_flags[i]);
  }

  Index intervalcount;
  intervals.getIntervalCount(intervalcount);
  for (Index i=0;i<intervalcount;i++) {






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

}

SplitChoice SmallerBranching::branchInterval(Index intervalID, SmartPtr<IpoptApplication> app)
{

}

std::vector<SplitChoice> AbsoluteLargerBranching::branchSensitivityMatrix(SmartPtr<IpoptApplication> app)
{

}

SplitChoice AbsoluteLargerBranching::branchInterval(Index intervalID, SmartPtr<IpoptApplication> app)
{

}

Intervaluation* assignIntervaluationMethod(SmartPtr<OptionsList> options)
{
  std::string branchvalue;
  if (options->GetStringValue("branchvalue",branchvalue ,"")){
    printf("\branchvalue is set to %s\n",branchvalue.c_str());
    if (branchvalue =="bound")
      return new OneBoundIntervaluation();
     else if (branchvalue=="product")
       return new BothBoundIntervaluation();
  } else {
    printf("\nNo branchvalue given!\n");
  }
    return new OneBoundIntervaluation();

}

SplitChoice OneBoundIntervaluation::intervaluateInterval(Index intervalID,SmartPtr<IpoptApplication> app)
{

}

SplitChoice BothBoundIntervaluation::intervaluateInterval(Index intervalID,SmartPtr<IpoptApplication> app)
{

}

virtual SplitDecision SelectFirstControl::decideControlMethod(std::vector<SplitChoice> choices)
{
  SplitDecision retval;
  retval.intervalID = choices[0].intervalID;
  retval.parameterID = choices[0].intervalID;
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
  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(sens_matrix));

  const Index nrows = mv_sens->NRows();
  const Index ncols = mv_sens->NCols();

  // vector to store the indexes of sensitivity matrix rows containing control related data
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
  //  std::vector<Number> intervalwidths;
  SmartPtr<OptionsList> options = app->Options();
  //  IntervalWidthScaling* scalingmode = assignScalingMethod(options);
  //  scalingmode->scaleIntervalWidths(intervals,intervalwidths);
  BranchingCriterion* branchmode = assignBranchingMethod(options);
  //  branchmode->branchSensitivityMatrix();

  //  for (int i=0; i<intervalwidths.size();i++)
  // printf("\n skalierte Intervalbreite %d beträgt: %f\n",i,intervalwidths[i]);


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

  std::string branch_mode;
  if (options->GetStringValue("branchmode",branch_mode ,""))
    printf("\nbranch_mode is set to %s\n",branch_mode.c_str());
  else printf("\nNo branch_mode given!\n");

  std::string branchvalue;
  if (options->GetStringValue("branchvalue",branchvalue ,""))
    printf("\nbranchvalue is set to %s\n",branchvalue.c_str());
  else printf("\nNo branchvalue given!\n");

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
/*
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
}*/

Number Abs(Number value)
{
  Number retval =0;
  if (value<0)
    retval = (-value);
  else
    retval = value;
  return retval;
}

