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
#include "IpExpansionMatrix.hpp"

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
  IntervalInfoSet();
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
  //IntervalInfoSet();
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

struct SplitApproximation
{
  Index intervalID;
  Index parameterID;
  std::vector<Number> approx_sensitivities;
  // propably more to come (indices)
};

class BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<IpoptApplication> app) const = 0;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const = 0;
private:
  std::vector<Index> ctrl_rows_;
  SmartPtr<Matrix> sens_matrix_;
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

class SplitAlgorithm
{
public:
  virtual SplitDecision applySplitAlgorithm(SmartPtr<IpoptApplication> app) = 0;
};

class SplitWRTControlSensitivities : public SplitAlgorithm
{
  SplitDecision applySplitAlgorithm(SmartPtr<IpoptApplication> app);
};

class LinearizeKKTWithMINRES : public SplitAlgorithm
{
public:
  LinearizeKKTWithMINRES(SmartPtr<IpoptApplication> app);
  Index computeRHSDim() const;
  SplitDecision applySplitAlgorithm(SmartPtr<IpoptApplication> app);
  SplitApproximation applyAlgorithmOnInterval(SmartPtr<IpoptApplication> app, const Index& interval);
  bool splitIntervalIndices(const std::vector<Index>& intervalIDs,std::vector<Index>& static_indices,std::vector<Index>& shift_indices,const Index& interval);
  bool assignIntAndParaDepAttr(const Index& interval,const Index& column);
  SmartPtr<const DenseVector> extractColumn(SmartPtr<const Matrix> original,const Index& column) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<const DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<const DenseVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<const DenseVector> original,const std::vector<Index>& indices) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<DenseVector> original,const std::vector<Index>& indices) const;
  Number computeRMultS() const;
  Number computeSMultS() const;
  /*  SmartPtr<const DenseVector> computeR() const;*/ // not neccessary for x_0=0 and 0 steps.
  SmartPtr<const DenseVector> computeS() const;
  SplitDecision decideInterval(const std::vector<SplitApproximation>& approximates) const;
  /*  std::vector<SplitChoice>getChoicesFromApprox(const std::vector<SplitApproximation>& approx) const;
  const Number* getExpandedRHSValues() const;
*/
private:
  //original Ipopt data
  SmartPtr<const DenseVector> x_;
  SmartPtr<const DenseVector> p_;
  SmartPtr<const DenseVector> y_c_;
  SmartPtr<const DenseVector> y_d_;
  SmartPtr<OptionsList> options_;
  SmartPtr<const Matrix> rhs_h_;
  SmartPtr<const Matrix> rhs_c_;
  SmartPtr<const Matrix> rhs_d_;
  SmartPtr<const SymMatrix> lhs_h_;
  SmartPtr<const Matrix> lhs_c_;
  SmartPtr<const Matrix> lhs_d_;

  // extracted data - not interval specific

  IntervalInfoSet intervals_;
  std::vector<Index> u_indices_;
  std::vector<Index> x_intervalIDs_;
  std::vector<Index> p_intervalIDs_;
  std::vector<Index> c_intervalIDs_;
  std::vector<Index> d_intervalIDs_;
  Index n_i_;
  Index n_u_;
  Index n_p_;
  Index n_x_;
  Index n_c_;
  Index n_d_;
  Index rhs_dim_;

  // extracted data - interval specific
  SmartPtr<const DenseVector> rhs_static_h_;
  SmartPtr<const DenseVector> rhs_static_c_;
  SmartPtr<const DenseVector> rhs_static_d_;
  SmartPtr<const DenseVector> rhs_u_;
  SmartPtr<const DenseVector> rhs_shift_h_;
  SmartPtr<const DenseVector> rhs_shift_c_;
  SmartPtr<const DenseVector> rhs_shift_d_;
  std::vector<Index> static_x_indices_;
  std::vector<Index> shift_x_indices_;
  std::vector<Index> static_c_indices_;
  std::vector<Index> shift_c_indices_;
  std::vector<Index> static_d_indices_;
  std::vector<Index> shift_d_indices_;

  // manipulated or constructed data - interval specific
  SmartPtr<const DenseVector> rhs_i_;


  /*  // testwise
  SmartPtr<const DenseVector> rhs_h_i;
  SmartPtr<const DenseVector> rhs_c_i;
  SmartPtr<const DenseVector> rhs_d_i;
  */

};

SplitAlgorithm* assignSplitAlgorithm(SmartPtr<IpoptApplication> app);

//ParameterShift* assignShiftMethod(SmartPtr<OptionsList> options);

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
  // printf("\nAufruf des IntervalInfoSet::IntervalInfoSet(SmartPtr<const DenseVector> parameters) Konstruktors.");
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
  //  printf("\nLargerBranching::branchSensitivityMatrix: ctrl_rows.size() is: %d\n",ctrl_rows.size());
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
  //  printf("\nSmallerBranching::branchSensitivityMatrix: ctrl_rows.size() is: %d\n",ctrl_rows.size());
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

  //  printf("\nAbsoluteLargerBranching::branchSensitivityMatrix: ctrl_rows.size() is: %d\n",ctrl_rows.size());

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


SplitDecision SplitWRTControlSensitivities::applySplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  SmartPtr<OptionsList> options = app->Options();
  BranchingCriterion* branchmode = assignBranchingMethod(options);
  std::vector<SplitChoice> splitchoices = branchmode->branchSensitivityMatrix(app);

  ControlSelector* pickfirst = assignControlMethod(options);
  return pickfirst->decideSplitControl(splitchoices);
}

LinearizeKKTWithMINRES::LinearizeKKTWithMINRES(SmartPtr<IpoptApplication> app)
{
  x_intervalIDs_.clear();
  p_intervalIDs_.clear();
  c_intervalIDs_.clear();
  d_intervalIDs_.clear();

  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  //assign local original Ipopt data
  x_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->x()));
  p_ = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  y_c_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->y_c()));
  y_d_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->y_d()));
  options_ = app->Options();
  rhs_h_ = orig_nlp->h_p(*x_, 1.0, *y_c_, *y_d_);
  rhs_c_  = orig_nlp->jac_c_p(*x_);
  rhs_d_  = orig_nlp->jac_d_p(*x_);
  lhs_h_ = orig_nlp->h(*x_, 1.0, *y_c_, *y_d_);
  lhs_c_ = orig_nlp->jac_c(*x_);
  lhs_d_ = orig_nlp->jac_d(*x_);

  SmartPtr<const DenseVectorSpace> x_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_->OwnerSpace()));
  SmartPtr<const DenseVectorSpace> c_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(y_c_->OwnerSpace()));
  SmartPtr<const DenseVectorSpace> d_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(y_d_->OwnerSpace()));

  //assign local extracted data - non interval specific
  //  printf("\nES WIRD DER RICHTIGE KONSTRUKTOR VERWENDET!");
  IntervalInfoSet intervals_ = IntervalInfoSet(p_);
  //  printf("\nES WIRD DER RICHTIGE KONSTRUKTOR VERWENDET!");
  std::vector<Index> x_intervalIDs = x_space->GetIntegerMetaData("intervalID");
  std::vector<Index> u_indices;
  std::vector<Index> x_indices;
  if (splitIntervalIndices(x_intervalIDs,x_indices,u_indices,0)) {
    assert(u_indices.size());
    assert(x_indices.size());
    u_indices_ = u_indices;
    x_intervalIDs_=x_intervalIDs;



    }

  p_intervalIDs_ = intervals_.getIntervalIDs();
  if (c_space->HasIntegerMetaData("intervalID")){
    c_intervalIDs_ = c_space->GetIntegerMetaData("intervalID");
    // for (int i=0;i<c_intervalIDs_.size();i++) {
    //   printf("LinearizeKKTWithMINRES::LinearizeKKTWithMINRES(): c_intervalIDs_[%d]=%d\n",i,c_intervalIDs_[i]);
    // }
  } else
    printf("LinearizeKKTWithMINRES::LinearizeKKTWithMINRES(): Error - no intervalIDs for c_space!\n");
  if (d_space->HasIntegerMetaData("intervalID")) {
    d_intervalIDs_ = d_space->GetIntegerMetaData("intervalID");
    // for (int i=0;i<d_intervalIDs_.size();i++) {
    //   printf("LinearizeKKTWithMINRES::LinearizeKKTWithMINRES(): d_intervalIDs_[%d]=%d\n",i,d_intervalIDs_[i]);
    // }
  } else
    printf("LinearizeKKTWithMINRES::LinearizeKKTWithMINRES(): Error - no intervalIDs for d_space!\n");

  n_i_ = intervals_.getIntervalCount();
  n_u_ = u_indices_.size();
  n_p_ = intervals_.Size();
  n_x_ = x_intervalIDs_.size()-n_u_;
  n_c_ = c_intervalIDs_.size();
  n_d_ = d_intervalIDs_.size();

  rhs_dim_ = computeRHSDim();

  // interval specific extracted data cannot be initialized in constructor!!

}

/* apply MINRES-approximation of a split to all intervals and decide for the smartest split*/
SplitDecision LinearizeKKTWithMINRES::applySplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  std::vector<SplitApproximation> approximates(n_i_);
  printf("\nLinearizeKKTWithMINRES::applySplitAlgorithm(): n_i_ = %d",n_i_);
  for (int i=0;i<n_i_;i++) {

      // get MINRES approximated split results wrt each single interval
      // i+1: in case intervalIDs start with 1 (which they do)
      approximates[i] = applyAlgorithmOnInterval(app,i+1);
  }

  SplitDecision retval;

  // chose the split with best results
  retval = this->decideInterval(approximates);


  ///////////////////////only for the sake of a working python/ampl interface////
  SmartPtr<OptionsList> options = app->Options();
  BranchingCriterion* branchmode = assignBranchingMethod(options);
  std::vector<SplitChoice> splitchoices = branchmode->branchSensitivityMatrix(app);
  ControlSelector* pickfirst = assignControlMethod(options);
  retval = pickfirst->decideSplitControl(splitchoices);

  return retval;
}

Index LinearizeKKTWithMINRES::computeRHSDim() const
{
  Index retval=0;
  Index tmp_ID=0;
  Index tmp_i=0;

  // determine a feasible intervalID
  for (int i=0;i<x_intervalIDs_.size();i++) {
    if (x_intervalIDs_[i]) {
      tmp_ID=x_intervalIDs_[i];
      tmp_i=i;
      break;
    }
  }
  // add shift_x size
  for (int i=tmp_i;i<x_intervalIDs_.size();i++) {
    if (x_intervalIDs_[i]==tmp_ID)
      retval++;
  }
  // add shift_c size
  for (int i=0;i<c_intervalIDs_.size();i++) {
    if (c_intervalIDs_[i]==tmp_ID)
      retval++;
  }
  // add shift_d size
  for (int i=tmp_i;i<d_intervalIDs_.size();i++) {
    if (d_intervalIDs_[i]==tmp_ID)
      retval++;
  }
  //  printf("Index LinearizeKKTWithMINRES::computeRHSDim(): retval = %d,   n_x_ = %d, n_c_ = %d, n_d_ = %d,   n_u_ = %d ",retval, n_x_, n_c_, n_d_, n_u_);
  retval += n_x_ +n_c_ +n_d_+n_u_;

  return retval;
}

bool LinearizeKKTWithMINRES::assignIntAndParaDepAttr(const Index& interval,const Index& column)
{
  bool retval = 1;

    // split interval indices into to be shifted interval entries and remainder
  if (!splitIntervalIndices(x_intervalIDs_,static_x_indices_,shift_x_indices_,interval)) {
    printf("\nLinearizeKKTWithMINRES::assignIntAndParaDepAttr(): ERROR: unable to split x_intervalIDs_");
    retval = 0;
  }
  /*  for (int i=0;i<static_x_indices_.size();i++)
      printf("\nstatic_x_indices_[%d] = %d",i,static_x_indices_[i]);
  printf("\n\n");

  */

  if (!splitIntervalIndices(c_intervalIDs_,static_c_indices_,shift_c_indices_,interval)) {
    printf("\nLinearizeKKTWithMINRES::assignIntAndParaDepAttr(): ERROR: unable to split c_intervalIDs_");
    retval = 0;
  }
  if (!splitIntervalIndices(d_intervalIDs_,static_d_indices_,shift_d_indices_,interval)) {
    printf("\nLinearizeKKTWithMINRES::assignIntAndParaDepAttr(): ERROR: unable to split d_intervalIDs_");
    retval = 0;
  }
  // get different rhs entries for this column
  SmartPtr<const DenseVector> rhs_h_i;
  if (GetRawPtr(rhs_h_))
    rhs_h_i = extractColumn(rhs_h_,column);
  else
    printf("\nLinearizeKKTWithMINRES::assignIntAndParaDepAttr(): ERROR: rhs_h_ is NULL");
  SmartPtr<const DenseVector> rhs_c_i;
  if (GetRawPtr(rhs_c_))
    rhs_c_i = extractColumn(rhs_c_,column);
  else
    printf("\nLinearizeKKTWithMINRES::assignIntAndParaDepAttr(): ERROR: rhs_c_ is NULL");
  SmartPtr<const DenseVector> rhs_d_i;
  if (GetRawPtr(rhs_d_))
    rhs_d_i = extractColumn(rhs_d_,column);
  else
    printf("\nLinearizeKKTWithMINRES::assignIntAndParaDepAttr(): ERROR: rhs_d_ is NULL");

  // assign rhs subvectors
  //  printf("\n---------------------------------1");
  rhs_static_h_ = shrinkVector(rhs_h_i,static_x_indices_);
  //printf("\n---------------------------------2");
  rhs_static_c_ = shrinkVector(rhs_c_i,static_c_indices_);
//  printf("\n---------------------------------3");
  rhs_static_d_ = shrinkVector(rhs_d_i,static_d_indices_);
//  printf("\n---------------------------------4");
  rhs_shift_h_ = shrinkVector(rhs_h_i,shift_x_indices_);
//  printf("\n---------------------------------5");
  rhs_shift_c_ = shrinkVector(rhs_c_i,shift_c_indices_);
//  printf("\n---------------------------------6");
  rhs_shift_d_ = shrinkVector(rhs_d_i,shift_d_indices_);
//  printf("\n---------------------------------7");
  rhs_u_ = shrinkVector(rhs_h_i,u_indices_);

  // expand all vector parts to final dimension
  Index abs_pos = 0;
  SmartPtr<const DenseVector> rhs_x_static_part = expandVector(rhs_static_h_,rhs_dim_,abs_pos);
  abs_pos += static_x_indices_.size();
  SmartPtr<const DenseVector> rhs_x_shift_upart = expandVector(rhs_shift_h_,rhs_dim_,abs_pos);
  abs_pos += shift_x_indices_.size();
  SmartPtr<const DenseVector> rhs_u_part = expandVector(rhs_u_,rhs_dim_,abs_pos);
  abs_pos += n_u_;
  SmartPtr<const DenseVector> rhs_c_static_part = expandVector(rhs_static_c_,rhs_dim_,abs_pos);
  abs_pos += static_c_indices_.size();
  SmartPtr<const DenseVector> rhs_c_shift_upart = expandVector(rhs_shift_c_,rhs_dim_,abs_pos);
  abs_pos += shift_c_indices_.size();
  SmartPtr<const DenseVector> rhs_d_static_part = expandVector(rhs_static_d_,rhs_dim_,abs_pos);
  abs_pos += static_d_indices_.size();
  SmartPtr<const DenseVector> rhs_d_shift_upart = expandVector(rhs_shift_d_,rhs_dim_,abs_pos);
  abs_pos += shift_d_indices_.size();
  SmartPtr<const DenseVector> rhs_x_shift_lpart = expandVector(rhs_shift_h_,rhs_dim_,abs_pos);
  abs_pos += shift_x_indices_.size();
  SmartPtr<const DenseVector> rhs_c_shift_lpart = expandVector(rhs_shift_c_,rhs_dim_,abs_pos);
  abs_pos += shift_c_indices_.size();
  SmartPtr<const DenseVector> rhs_d_shift_lpart = expandVector(rhs_shift_d_,rhs_dim_,abs_pos);
  abs_pos += shift_d_indices_.size();

  // add all vectors
  SmartPtr<DenseVectorSpace> rhs_i_space = new DenseVectorSpace(rhs_dim_);
  SmartPtr<DenseVector> rhs_i = dynamic_cast<const DenseVector*>(rhs_i_space->MakeNewDenseVector());
  rhs_i->AddTwoVectors(1.0,*rhs_x_static_part,1.0,*rhs_c_static_part,0.0);
  rhs_i->AddTwoVectors(1.0,*rhs_d_static_part,1.0,*rhs_x_shift_upart,1.0);
  rhs_i->AddTwoVectors(1.0,*rhs_c_shift_upart,1.0,*rhs_d_shift_upart,1.0);
  rhs_i->AddTwoVectors(1.0,*rhs_x_shift_lpart,1.0,*rhs_u_part,1.0);
  rhs_i->AddTwoVectors(1.0,*rhs_c_shift_lpart,1.0,*rhs_d_shift_lpart,1.0);
  rhs_i_ = dynamic_cast<const DenseVector*>(GetRawPtr(rhs_i));

  return retval;
}

// get a specific column of a given const Matrix* as a const Densevector*
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::extractColumn(SmartPtr<const Matrix> original,const Index& column) const
{
  const Index n_cols = original->NCols();
  const Index n_rows = original->NRows();
  // ncols+1 because starting to count at 1...
  assert(column<n_cols+1);

  //setup unit vector with sufficient dimension
  SmartPtr<DenseVectorSpace> unit_space = new DenseVectorSpace(n_cols);
  SmartPtr<DenseVector> unit = dynamic_cast<const DenseVector*>(unit_space->MakeNewDenseVector());
  Number* unit_values = new Number[n_cols];
  for (int i=0;i<n_cols;++i)
    unit_values[i]=0.0;
  unit_values[column]=1.0;
  unit->SetValues(unit_values);

  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(n_rows);
  SmartPtr<DenseVector> retval = dynamic_cast<const DenseVector*>(retval_space->MakeNewDenseVector());
  original->MultVector(1.0,*unit,0.0,*retval);

  return dynamic_cast<const DenseVector*>(GetRawPtr(retval));
}

/* expand Vector original from smaller (old) dim to large_dim, with the vector indices listing at which indices in the new vector the old values are to be found */
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::expandVector(SmartPtr<const DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const
{
  const Index small_dim = indices.size();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d",small_dim,large_dim);
  assert(large_dim>=small_dim);

  if (large_dim==small_dim) {
//    printf("\nLinearizeKKTWithMINRES::expandVector(): WARNING: called to expand a vector to original size.");
    return original;
  } else {
//    printf("\n");
//    for (int j=0;j<indices.size();j++)
//      printf("\nindices[%d] = %d",j,indices[j]);
    Index* exppos = new Index[small_dim];
    for (int i=0;i<small_dim;i++)
      exppos[i]=indices[i];
//    printf("\n");
//    for (int i=0;i<small_dim;i++)
//      printf("\nexppos[%d] = %d",i,exppos[i]);

    SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
    SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
    SmartPtr<DenseVectorSpace> expanded_space = new DenseVectorSpace(large_dim);
    SmartPtr<Vector> new_vec = expanded_space->MakeNew();

    em->MultVector(1.0,*original,0.0,*new_vec);

    return dynamic_cast<const DenseVector*>(GetRawPtr(new_vec));
  }
  printf("\nLinearizeKKTWithMINRES::expandVector(): Unknown ERROR.");
  return NULL;
}

/*expand Vector original from smaller (old) dim to large_dim, inserting the original values as a dense block at start_idx in the new vector*/
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::expandVector(SmartPtr<const DenseVector> original,const Index& large_dim, const Index& start_idx) const
{
  const Index small_dim = original->Dim();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d   start_idx = %d",small_dim,large_dim,start_idx);
  assert(large_dim>=small_dim);

  if (large_dim==small_dim) {
    printf("\nLinearizeKKTWithMINRES::expandVector(): WARNING: called to expand a vector to original size.");
    return original;
  } else {
    Index* exppos = new Index[small_dim];
    for (int i=0;i<small_dim;i++)
      exppos[i]=i+start_idx;
//    printf("\n");
//    for (int i=0;i<small_dim;i++)
//      printf("\nexppos[%d] = %d",i,exppos[i]);

    SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
    SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
    SmartPtr<DenseVectorSpace> expanded_space = new DenseVectorSpace(large_dim);
    SmartPtr<Vector> new_vec = expanded_space->MakeNew();

    em->MultVector(1.0,*original,0.0,*new_vec);

    return dynamic_cast<const DenseVector*>(GetRawPtr(new_vec));
  }
  printf("\nLinearizeKKTWithMINRES::expandVector(): Unknown ERROR.");
  return NULL;
}

/* shrink Vector original from larger (old) dim to smaller dim, with the vector indices listing which indices from the old vector to keep*/
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::shrinkVector(SmartPtr<const DenseVector> original, const std::vector<Index>& indices) const
{
  const Index large_dim = original->Dim();
  const Index small_dim = indices.size();
//  printf("\nshrinkVector(): large_dim = %d   small_dim = %d",large_dim,small_dim);
  assert(large_dim>=small_dim);
  if (large_dim==small_dim) {
    printf("\nLinearizeKKTWithMINRES::shrinkVector(const): WARNING: called to shrink a vector to original size.");
    return original;
  } else {
//    printf("\n");
//    for (int j=0;j<indices.size();j++)
//      printf("\nindices[%d] = %d",j,indices[j]);
    Index* exppos = new Index[small_dim];
    for (int i=0;i<small_dim;i++)
      exppos[i]=indices[i];
//    printf("\n");
//    for (int i=0;i<small_dim;i++)
//      printf("\nexppos[%d] = %d",i,exppos[i]);

    SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
    SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
    SmartPtr<DenseVectorSpace> expanded_space = new DenseVectorSpace(small_dim);
    SmartPtr<Vector> new_vec = expanded_space->MakeNew();

    em->TransMultVector(1.0,*original,0.0,*new_vec);

    return dynamic_cast<const DenseVector*>(GetRawPtr(new_vec));
  }
  printf("\nLinearizeKKTWithMINRES::shrinkVector(const): Unknown ERROR.");
  return NULL;
}

/* same as shrinkVector, just with nonconst DenseVector entering*/
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::shrinkVector(SmartPtr<DenseVector> original, const std::vector<Index>& indices) const
{
  SmartPtr<const DenseVector> c_original = dynamic_cast<const DenseVector*>(GetRawPtr(original));

  return shrinkVector(c_original,indices);
}

Number LinearizeKKTWithMINRES::computeRMultS() const
{ /*
  Number retval;
  SmartPtr<const DenseVector> S = computeS();
  retval = rhs_i_->Dot(*S);

  return retval; */
}

Number LinearizeKKTWithMINRES::computeSMultS() const
{ /*
  Number retval;
  Number squareval;

  SmartPtr<const DenseVector> S = computeS();

  retval = S->Dot(*S);
  squareval = (S->Nrm2());
  squareval = squareval*squareval;
  if (retval==squareval)
    printf("\nLinearizeKKTWithMINRES::computeSMultS(): Both ways of computing SMultS work fine.");
    return retval; */
}

/* SmartPtr<const DenseVector> LinearizeKKTWithMINRES::computeR() const
{
// not needed for x_0 = 0 and 0 steps.
} */

SmartPtr<const DenseVector> LinearizeKKTWithMINRES::computeS() const
{
  // static non-control parts of S equal static rhs part
  SmartPtr<const DenseVector> retval_st_h = rhs_static_h_;
  SmartPtr<const DenseVector> retval_st_c = rhs_static_c_;
  SmartPtr<const DenseVector> retval_st_d = rhs_static_d_;

  // setup u part
  SmartPtr<DenseVector> retval_u = dynamic_cast<DenseVector*>(rhs_u_->MakeNew());

  // get (W_i,shift times y_shift) with i= (shift, u)
  //  SmartPtr<DenseVectorSpace> w_extr_space = new DenseVectorSpace(lhs_h_->NCols());
  SmartPtr<const DenseVector> y_shift_large = expandVector(rhs_shift_h_,shift_x_indices_,lhs_h_->NCols());
  SmartPtr<DenseVector> w_extractor = dynamic_cast<DenseVector*>(y_shift_large->MakeNew());
  lhs_h_->MultVector(1.0,*y_shift_large,0.0,*w_extractor);


  // map result back to get u and shift part seperately for the appropriate parts of S
  //  SmartPtr<const DenseVector> u_part_wy = shrinkVector(dynamic_cast);


  // upper shifted parts of S equal shift rhs part
  SmartPtr<const DenseVector> retval_ush_h = rhs_shift_h_;
  SmartPtr<const DenseVector> retval_ush_c =  rhs_shift_c_;
  SmartPtr<const DenseVector> retval_ush_d =  rhs_shift_d_;

  // setup lower shifted h part assigned to x
  SmartPtr<DenseVector> retval_lsh_h = dynamic_cast<DenseVector*>(rhs_shift_h_->MakeNew());





  // setup lower shifted c part
  SmartPtr<DenseVector> retval_lsh_c = dynamic_cast<DenseVector*>(rhs_shift_c_->MakeNew());
  // setup lower shifted d part
  SmartPtr<DenseVector> retval_lsh_d = dynamic_cast<DenseVector*>(rhs_shift_d_->MakeNew());





  return NULL;
}

SplitDecision LinearizeKKTWithMINRES::decideInterval(const std::vector<SplitApproximation>& approximates) const
{
  /* waiting to be implemeted*/
}

/*
SmartPtr<Vector> LinearizeKKTWithMINRES::getAshiftMultYshift(const Index& column,SmartPtr<IpoptApplication> app,const Index & interval,const RHSVectors& shifted_rhs) const
{
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  // get raw material for Ashift
  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<const Matrix> jac_c = orig_nlp->jac_c(*x);
  printf("\ngetAshiftMultYshift: jac_c->NCols()= %d, jac_c->NRows()= %d\n",jac_c->NCols(),jac_c->NRows());

  // set up ExpansionMatrix mapping upwards from rhs
  SmartPtr<const VectorSpace> x_space = x->OwnerSpace();
  const std::vector<Index> intervalIDs = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_space))->GetIntegerMetaData("intervalID");
  Index small_dim = shifted_rhs.shift_h_part->Dim();
  Index large_dim = x->Dim();
  Index* upw_exppos = new Index[small_dim];
  Index curr_pos = 0;
  for (int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i]==interval) {
      upw_exppos[curr_pos]=i;
      curr_pos++;
    }
  }
  SmartPtr<ExpansionMatrixSpace> upw_em_space = new ExpansionMatrixSpace(large_dim,small_dim,upw_exppos);
  SmartPtr<ExpansionMatrix> upw_em = upw_em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> upw_cache_space = new DenseVectorSpace(large_dim);
  SmartPtr<Vector> upw_cache = upw_cache_space->MakeNew();

  // set up ExpansionMatrix mapping downwards to return value
  small_dim = shifted_rhs.shift_c_part->Dim();
  large_dim = jac_c->NRows();
  Index* dnw_exppos = new Index[small_dim];
  for (int i=0;i<small_dim;i++)
    dnw_exppos[i] = i+(interval-1)*small_dim;
  SmartPtr<ExpansionMatrixSpace> dnw_em_space = new ExpansionMatrixSpace(large_dim,small_dim,dnw_exppos);
  SmartPtr<ExpansionMatrix> dnw_em = dnw_em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> dnw_cache_space = new DenseVectorSpace(large_dim);
  SmartPtr<Vector> dnw_cache = dnw_cache_space->MakeNew();

  // transform Yshift back and forth
  SmartPtr<Vector> retval = shifted_rhs.shift_c_part->MakeNew();
  upw_em->MultVector(1.0,*shifted_rhs.shift_h_part,0.0,*upw_cache);
  jac_c->MultVector(1.0,*upw_cache,0.0,*dnw_cache);
  //  printf("\ngetWuschiftMultYshift(): final extraction dimensions: retval->Dim() = %d, cache->Dim() = %d  em->NCols() = %d, em->NRows() = %d\n",retval->Dim(),cache->Dim(),em->NCols(),em->NRows());
  dnw_em->TransMultVector(1.0,*dnw_cache,0.0,*retval);

  return retval;
}

SmartPtr<Vector> LinearizeKKTWithMINRES::getBshiftMultYshift(const Index& column,SmartPtr<IpoptApplication> app,const Index & interval,const RHSVectors& shifted_rhs) const
{
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  // get raw material for Bshift
  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<const Matrix> jac_d = orig_nlp->jac_d(*x);
  printf("\ngetBshiftMultYshift: jac_d->NCols()= %d, jac_d->NRows()= %d\n",jac_d->NCols(),jac_d->NRows());

  // set up ExpansionMatrix mapping upwards from rhs
  SmartPtr<const VectorSpace> x_space = x->OwnerSpace();
  const std::vector<Index> intervalIDs = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_space))->GetIntegerMetaData("intervalID");
  Index small_dim = shifted_rhs.shift_h_part->Dim();
  Index large_dim = x->Dim();
  Index* upw_exppos = new Index[small_dim];
  Index curr_pos = 0;
  for (int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i]==interval) {
      upw_exppos[curr_pos]=i;
      curr_pos++;
    }
  }
  SmartPtr<ExpansionMatrixSpace> upw_em_space = new ExpansionMatrixSpace(large_dim,small_dim,upw_exppos);
  SmartPtr<ExpansionMatrix> upw_em = upw_em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> upw_cache_space = new DenseVectorSpace(large_dim);
  SmartPtr<Vector> upw_cache = upw_cache_space->MakeNew();

  // set up ExpansionMatrix mapping downwards to return value
  small_dim = shifted_rhs.shift_d_part->Dim();
  large_dim = jac_d->NRows();
  Index* dnw_exppos = new Index[small_dim];
  for (int i=0;i<small_dim;i++)
    dnw_exppos[i] = i+(interval-1)*small_dim;
  SmartPtr<ExpansionMatrixSpace> dnw_em_space = new ExpansionMatrixSpace(large_dim,small_dim,dnw_exppos);
  SmartPtr<ExpansionMatrix> dnw_em = dnw_em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> dnw_cache_space = new DenseVectorSpace(large_dim);
  SmartPtr<Vector> dnw_cache = dnw_cache_space->MakeNew();

  // transform Yshift back and forth
  SmartPtr<Vector> retval = shifted_rhs.shift_d_part->MakeNew();
  upw_em->MultVector(1.0,*shifted_rhs.shift_h_part,0.0,*upw_cache);
  jac_d->MultVector(1.0,*upw_cache,0.0,*dnw_cache);
  //  printf("\ngetWuschiftMultYshift(): final extraction dimensions: retval->Dim() = %d, cache->Dim() = %d  em->NCols() = %d, em->NRows() = %d\n",retval->Dim(),cache->Dim(),em->NCols(),em->NRows());
  dnw_em->TransMultVector(1.0,*dnw_cache,0.0,*retval);

  return retval;
}

SmartPtr<Vector> LinearizeKKTWithMINRES::getWushiftMultYshift(const Index& column,SmartPtr<IpoptApplication> app,const Index & interval,const RHSVectors& shifted_rhs) const
{
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  // get raw material for Wshift
  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<const Vector> y_c = app->IpoptDataObject()->curr()->y_c();
  SmartPtr<const Vector> y_d = app->IpoptDataObject()->curr()->y_d();
  SmartPtr<const SymMatrix> h = orig_nlp->h(*x, 1.0, *y_c, *y_d);
  printf("\ngetWushiftMultYshift: h->NCols()= %d, h->NRows()= %d\n",h->NCols(),h->NRows());

  // set up ExpansionMatrix mapping upwards from rhs
  SmartPtr<const VectorSpace> x_space = x->OwnerSpace();
  const std::vector<Index> intervalIDs = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_space))->GetIntegerMetaData("intervalID");
  Index small_dim = shifted_rhs.shift_h_part->Dim();
  Index large_dim = x->Dim();
  Index* upw_exppos = new Index[small_dim];
  Index curr_pos = 0;
  for (int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i]==interval) {
      upw_exppos[curr_pos]=i;
      curr_pos++;
    }
  }
  SmartPtr<ExpansionMatrixSpace> upw_em_space = new ExpansionMatrixSpace(large_dim,small_dim,upw_exppos);
  SmartPtr<ExpansionMatrix> upw_em = upw_em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> upw_cache_space = new DenseVectorSpace(large_dim);
  SmartPtr<Vector> upw_cache = upw_cache_space->MakeNew();

  // set up ExpansionMatrix mapping downwards to return value
  small_dim = shifted_rhs.u_part->Dim();
  Index* dnw_exppos = new Index[small_dim];
  curr_pos = 0;
  for (int i=0;i<intervalIDs.size();i++) {
    if (!intervalIDs[i]) {
      dnw_exppos[curr_pos]=i;
      curr_pos++;
    }
  }
  SmartPtr<ExpansionMatrixSpace> dnw_em_space = new ExpansionMatrixSpace(large_dim,small_dim,dnw_exppos);
  SmartPtr<ExpansionMatrix> dnw_em = dnw_em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> dnw_cache_space = new DenseVectorSpace(large_dim);
  SmartPtr<Vector> dnw_cache = dnw_cache_space->MakeNew();

  // transform Yshift back and forth
  SmartPtr<Vector> retval = shifted_rhs.u_part->MakeNew();
  upw_em->MultVector(1.0,*shifted_rhs.shift_h_part,0.0,*upw_cache);
  h->MultVector(1.0,*upw_cache,0.0,*dnw_cache);
  //  printf("\ngetWuschiftMultYshift(): final extraction dimensions: retval->Dim() = %d, cache->Dim() = %d  em->NCols() = %d, em->NRows() = %d\n",retval->Dim(),cache->Dim(),em->NCols(),em->NRows());
  dnw_em->TransMultVector(1.0,*dnw_cache,0.0,*retval);

  return retval;
}

SmartPtr<Vector> LinearizeKKTWithMINRES::getWshiftMultYshift(const Index& column,SmartPtr<IpoptApplication> app,const Index & interval,const RHSVectors& shifted_rhs) const
{
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  // get raw material for Wshift
  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<const Vector> y_c = app->IpoptDataObject()->curr()->y_c();
  SmartPtr<const Vector> y_d = app->IpoptDataObject()->curr()->y_d();
  SmartPtr<const SymMatrix> h = orig_nlp->h(*x, 1.0, *y_c, *y_d);
  printf("\ngetWshiftMultYshift: h->NCols()= %d, h->NRows()= %d\n",h->NCols(),h->NRows());


  // setup ExpansionMatrix
  SmartPtr<const VectorSpace> x_space = x->OwnerSpace();
  const std::vector<Index> intervalIDs = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_space))->GetIntegerMetaData("intervalID");
  Index small_dim = shifted_rhs.shift_h_part->Dim();
  Index large_dim = h->NRows();
  Index* exppos = new Index[small_dim];
  Index curr_pos = 0;
  for (int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i]==interval) {
      exppos[curr_pos]=i;
      curr_pos++;
    }
  }

  SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
  SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> cache_space = new DenseVectorSpace(large_dim);
  SmartPtr<Vector> cache = cache_space->MakeNew();

  em->MultVector(1.0,*shifted_rhs.shift_h_part,0.0,*cache);
  h->MultVector(1.0,*cache,0.0,*cache);
  SmartPtr<Vector> retval = shifted_rhs.shift_h_part->MakeNew();
  printf("\ngetWschiftMultYshift(): final extraction dimensions: retval->Dim() = %d, cache->Dim() = %d  em->NCols() = %d, em->NRows() = %d\n",retval->Dim(),cache->Dim(),em->NCols(),em->NRows());
  em->TransMultVector(1.0,*cache,0.0,*retval);

  return retval;
}

SmartPtr<Vector> LinearizeKKTWithMINRES::getDMultYshift(const Index& column,SmartPtr<IpoptApplication> app,const Index & interval,const RHSVectors& shifted_rhs) const
{
  // setup ExpansionMatrix
  SmartPtr<Matrix> sense = getSensitivityMatrix(app);
  Index small_dim = shifted_rhs.shift_h_part->Dim();
  Index large_dim = sense->NRows();
  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<const VectorSpace> x_space = x->OwnerSpace();
  const std::vector<Index> intervalIDs = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_space))->GetIntegerMetaData("intervalID");
  Index* exppos = new Index[small_dim];
  Index curr_pos = 0;
  for (int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i]==interval) {
      exppos[curr_pos]=i;
      curr_pos++;
    }
  }

  SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
  SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> cache_space = new DenseVectorSpace(large_dim);
  SmartPtr<Vector> cache = cache_space->MakeNew();

  em->MultVector(1.0,*shifted_rhs.shift_h_part,0.0,*cache);
  SmartPtr<Vector> retval = shifted_rhs.u_part->MakeNew();
  em->TransMultVector(1.0,*cache,0.0,*retval);

  return retval;
}
*/
SplitApproximation LinearizeKKTWithMINRES::applyAlgorithmOnInterval(SmartPtr<IpoptApplication> app, const Index& interval)
{
  // assign private attributes depending on chosen interval
  printf("\nLinearizeKKTWithMINRES::applyAlgorithmOnInterval(): n_p_ = %d\n", n_p_);
  for (Index j=0;j<n_p_;j++) {
    if (assignIntAndParaDepAttr(interval,j)) {
      // compute the solution of the MINRES-Step
      const Number steplength = computeRMultS()*computeSMultS();
      /*      printf("\n\n");
      char buffer[16];
      sprintf(buffer,"rhs_i_[%d]",j);
      rhs_i_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, buffer);
      */
      SmartPtr<Vector> lhs = rhs_i_->MakeNewCopy();

      /*
      printf("\n\n");
  for (int i=0;i<static_x_indices_.size();i++)
    printf("\n static_x_indices_[%d] = %d",i,static_x_indices_[i]);
      printf("\n\n");
  for (int i=0;i<shift_x_indices_.size();i++)
    printf("\n shift_x_indices_[%d] = %d",i,shift_x_indices_[i]);
      printf("\n\n");
  rhs_h_i->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_h_i");
      printf("\n\n");
  for (int i=0;i<u_indices_.size();i++)
    printf("\n u_indices_[%d] = %d",i,u_indices_[i]);
      printf("\n\n");
  rhs_u_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_u_");
      printf("\n\n");
  for (int i=0;i<static_c_indices_.size();i++)
    printf("\n static_c_indices_[%d] = %d",i,static_c_indices_[i]);
      printf("\n\n");
  for (int i=0;i<shift_c_indices_.size();i++)
    printf("\n shift_c_indices_[%d] = %d",i,shift_c_indices_[i]);
      printf("\n\n");
  rhs_c_i->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_c_i");
      printf("\n\n");
  for (int i=0;i<static_d_indices_.size();i++)
    printf("\n static_d_indices_[%d] = %d",i,static_d_indices_[i]);
      printf("\n\n");
  for (int i=0;i<shift_d_indices_.size();i++)
    printf("\n shift_d_indices_[%d] = %d",i,shift_d_indices_[i]);
      printf("\n\n");
  rhs_d_i->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_d_i");
      */

  //  printf("\n n_u_= %d, n_x_ = %d, n_c_ = %d, n_d_ = %d", n_u_, n_x_, n_c_, n_d_);

  //  rhs_static_h_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_static_h");
  // rhs_static_c_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_static_c");
  // rhs_static_d_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_static_d");
  //  rhs_u_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_u");
  //  rhs_shift_h_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_shift_h");
  // rhs_shift_c_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_shift_c");
  // rhs_shift_d_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"rhs_shift_d");
        //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (steplength)
	lhs->Scal(1/steplength);
      else
	printf("\nLinearizeKKTWithMINRES::applyAlgorithmOnInterval(): ERROR: MINRES-Step is 0.");
    }
    else
      printf("\nLinearizeKKTWithMINRES::applyAlgorithmOnInterval(): ERROR: Unable to assign interval dependent attributes!");

    //SmartPtr<MultiVectorMatrixSpace> sense_space = new MultiVectorMatrixSpace();
    // add the resulting sense vector to matrix ....
    // calculate absolute sensitivity values
  }

  /* // Mapping test by checking for 0s in the right spots
  SmartPtr<DenseVector> retval_lsh_h = dynamic_cast<DenseVector*>(rhs_shift_h_->MakeNew());

    retval_lsh_h->Set(1.0);
  SmartPtr<const DenseVector> test_vec = expandVector(dynamic_cast<const DenseVector*>(GetRawPtr(retval_lsh_h)),static_x_indices_,lhs_h_->NCols());
  SmartPtr<DenseVector> res_vec = dynamic_cast<DenseVector*>(test_vec->MakeNew());
  printf("\n\n");
  lhs_h_->MultVector(1.0,*test_vec,0.0,*res_vec);
  test_vec = shrinkVector(dynamic_cast<const DenseVector*>(GetRawPtr(res_vec)),shift_x_indices_);
  test_vec->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"first 0s??");

  printf("\n\n");
  test_vec = expandVector(dynamic_cast<const DenseVector*>(GetRawPtr(retval_lsh_h)),shift_x_indices_,lhs_h_->NCols());
  lhs_h_->MultVector(1.0,*test_vec,0.0,*res_vec);
  test_vec = shrinkVector(dynamic_cast<const DenseVector*>(GetRawPtr(res_vec)),static_x_indices_);
  test_vec->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"second 0s??");


  retval_lsh_h = dynamic_cast<DenseVector*>(rhs_shift_c_->MakeNew());
  retval_lsh_h->Set(1.0);

  printf("\n\n");
  test_vec = expandVector(dynamic_cast<const DenseVector*>(GetRawPtr(retval_lsh_h)),shift_x_indices_,lhs_c_->NCols());
  lhs_c_->MultVector(1.0,*test_vec,0.0,*res_vec);
  test_vec = shrinkVector(dynamic_cast<const DenseVector*>(GetRawPtr(res_vec)),static_c_indices_);
  test_vec->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"third 0s??");

  printf("\n\n");
  test_vec = expandVector(dynamic_cast<const DenseVector*>(GetRawPtr(retval_lsh_h)),static_x_indices_,lhs_c_->NCols());
  lhs_c_->MultVector(1.0,*test_vec,0.0,*res_vec);
  test_vec = shrinkVector(dynamic_cast<const DenseVector*>(GetRawPtr(res_vec)),shift_c_indices_);
  test_vec->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"fourth 0s??");
*/
    //  lhs_h_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "hhhhhhhhh");

  SplitApproximation retval;

  return retval;
}

/** given an Index-type vector containing a list of intervalIDs, loads indices matching the given interval into the shift_indices part, the rest into static_indices **/
bool LinearizeKKTWithMINRES::splitIntervalIndices(const std::vector<Index>& intervalIDs,std::vector<Index>& static_indices,std::vector<Index>& shift_indices,const Index& interval)
{
  // make sure nothings in the way of push_backs
  static_indices.clear();
  shift_indices.clear();

  assert(intervalIDs.size());

  // cycle through intervalIDs and assign the indices
  for (int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i] == interval)
      shift_indices.push_back(i);
    else if (intervalIDs[i])
	  static_indices.push_back(i);
  }
/*   printf("\n\n");
   for (int i=0;i<static_indices.size();i++)
     printf("\nstatic_indices[%d] = %d",i,static_indices[i]);
   for (int i=0;i<shift_indices.size();i++)
     printf("\nshift_indices[%d] = %d",i,shift_indices[i]);
   printf("\n\n");
*/

  return 1;
}

SplitAlgorithm* assignSplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  return new LinearizeKKTWithMINRES(app);
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
  // Add the suffixes for interval-organization
  suffix_handler->AddAvailableSuffix("intervalID", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("intervalID", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);

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

  SplitAlgorithm* splitter = assignSplitAlgorithm(app);
  SplitDecision resulting_split = splitter->applySplitAlgorithm(app);

  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();

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

