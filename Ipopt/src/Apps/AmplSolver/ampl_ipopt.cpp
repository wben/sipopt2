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

struct VectorSet
{
  SmartPtr<IteratesVector> top;
  SmartPtr<DenseVector> x;
  SmartPtr<DenseVector> y_c;
  SmartPtr<DenseVector> y_d;
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

class SelectIndexControl : public ControlSelector
{
public:
  SelectIndexControl(const Index& index);
  SplitDecision decideSplitControl(const std::vector<SplitChoice>& choices) const;
private:
  Index index_;
};

class SelectNameControl : public ControlSelector
{
public:
  SelectNameControl(const std::string& name);
  SplitDecision decideSplitControl(const std::vector<SplitChoice>& choices) const;
private:
  std::string name_;
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
  SmartPtr<const DenseVector> expandVector(SmartPtr<const Vector> original,const std::vector<Index>& indices, const Index& large_dim) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<const DenseVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<const Vector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<const DenseVector> original,const std::vector<Index>& indices) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<DenseVector> original,const std::vector<Index>& indices) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<const Vector> original,const std::vector<Index>& indices) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<Vector> original,const std::vector<Index>& indices) const;

  Number computeRMultS(SmartPtr<IpoptApplication> app) const;
  Number computeSMultS(SmartPtr<IpoptApplication> app) const;
  /*  SmartPtr<const DenseVector> computeR() const;*/ // not neccessary for x_0=0 and 0 steps.
  SmartPtr<const DenseVector> computeS(SmartPtr<IpoptApplication> app) const;
  SmartPtr<const DenseVector> computeDMultYrhs(SmartPtr<IpoptApplication> app) const;
  SmartPtr<const DenseVector> undoConditionning(SmartPtr<IpoptApplication> app, SmartPtr<Vector> cond_sol) const;
  SmartPtr<const DenseVector> computeDeltaP(const IntervalInfoSet& intervals,const Index& interval) const;
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
  Index n_st_x_;
  Index n_sh_x_;
  Index n_st_c_;
  Index n_sh_c_;
  Index n_st_d_;
  Index n_sh_d_;

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

class LinKKTFaster : public SplitAlgorithm
{
public:
  LinKKTFaster(SmartPtr<IpoptApplication> app);
  Index computeRHSDim() const;
  SplitDecision applySplitAlgorithm(SmartPtr<IpoptApplication> app);
  SplitApproximation applyAlgorithmOnInterval(SmartPtr<IpoptApplication> app, const Index& interval);
  bool splitIntervalIndices(const std::vector<Index>& intervalIDs,std::vector<Index>& static_indices,std::vector<Index>& shift_indices,const Index& interval);
  bool assignIntAndParaDepAttr(const Index& interval,const Index& column);
  SmartPtr<const DenseVector> extractColumn(SmartPtr<const Matrix> original,const Index& column) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<const DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<const Vector> original,const std::vector<Index>& indices, const Index& large_dim) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<const DenseVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<const DenseVector> expandVector(SmartPtr<const Vector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<const DenseVector> original,const std::vector<Index>& indices) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<DenseVector> original,const std::vector<Index>& indices) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<const Vector> original,const std::vector<Index>& indices) const;
  SmartPtr<const DenseVector> shrinkVector(SmartPtr<Vector> original,const std::vector<Index>& indices) const;
  SmartPtr<DenseVector> expand(SmartPtr<DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const;
  SmartPtr<DenseVector> expand(SmartPtr<DenseVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<DenseVector> expand(SmartPtr<IteratesVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<DenseVector> shrink(SmartPtr<DenseVector> original,const std::vector<Index>& indices) const;
  SmartPtr<DenseVector> resetValuesAt(SmartPtr<DenseVector> target, std::vector<Index>indices) const;
  Number computeRMultS(SmartPtr<IpoptApplication> app,const VectorSet& target) const;
  Number computeSMultS(SmartPtr<IpoptApplication> app,const VectorSet& target) const;
  /*  SmartPtr<const DenseVector> computeR() const;*/ // not neccessary for x_0=0 and 0 steps.
  VectorSet computeS(SmartPtr<IpoptApplication> app,const VectorSet& target) const;
  VectorSet computeDMultYrhs(SmartPtr<IpoptApplication> app,const VectorSet& target) const;
  VectorSet undoConditionning(SmartPtr<IpoptApplication> app, const VectorSet& target) const;
  SmartPtr<const DenseVector> computeDeltaP(const IntervalInfoSet& intervals,const Index& interval) const;
  SplitDecision decideInterval(const std::vector<SplitApproximation>& approximates) const;
  VectorSet scaleVectorSet(const VectorSet& target,const Number& factor) const;
  SmartPtr<DenseVector> transformVectorSet(const VectorSet& target) const;
  void printVectorSet(SmartPtr<IpoptApplication>app,const VectorSet& target, const std::string& label) const;
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
  Index n_st_x_;
  Index n_sh_x_;
  Index n_st_c_;
  Index n_sh_c_;
  Index n_st_d_;
  Index n_sh_d_;

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
  SmartPtr<const DenseVector> x_i_;
  SmartPtr<const DenseVector> y_c_i_;
  SmartPtr<const DenseVector> y_d_i_;


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
  Index nints = getIntervalCount();
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
  if (options->GetStringValue("branchmode",branchmode ,"")) {
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

SelectIndexControl::SelectIndexControl(const Index& index)
{
  index_ = index;
}

SplitDecision SelectIndexControl::decideSplitControl(const std::vector<SplitChoice>& choices) const
{
  assert(choices.size()-index_>0);
  //simply returns index_th entry of SplitDecisions
  SplitDecision retval;
  retval.intervalID = choices[index_].intervalID;
  retval.parameterID = choices[index_].parameterID;
  return retval;
}

SelectNameControl::SelectNameControl(const std::string& name)
{
  name_ = name;
}

SplitDecision SelectNameControl::decideSplitControl(const std::vector<SplitChoice>& choices) const
{
  /*    //////////////////////////////////////////////////not implemented yet
  SplitDecision retval;
  retval.intervalID = choices[index_].intervalID;
  retval.parameterID = choices[index_].parameterID;
  return retval; */
}



ControlSelector* assignControlMethod(SmartPtr<OptionsList> options)
{
  // read from options what kind of controlselector to return
  ControlSelector* retval = new SelectIndexControl(0);
  Index index;
  std::string name;
  if (options->GetStringValue("ctrl_name",name,"")) {
    if (name.size()>0)
      retval =  new SelectNameControl(name);
  } else if (options->GetIntegerValue("ctrl_index",index,"")) {

    retval =  new SelectIndexControl(index);
  } else {
    printf("\nassignControlMethod(); The options don't assign a control for selection!\n");
  }
  return retval;
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
  intervals_ = IntervalInfoSet(p_);
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

  n_sh_x_ = int(n_x_/n_i_);
  n_st_x_ = n_x_ - n_sh_x_;
  n_sh_c_ = int(n_c_/n_i_);
  n_st_c_ = n_c_ - n_sh_c_;
  n_sh_d_ = int(n_d_/n_i_);
  n_st_d_ = n_d_ - n_sh_d_;

  rhs_dim_ = computeRHSDim();

  // interval specific extracted data cannot be initialized in constructor!!

}

/* apply MINRES-approximation of a split to all intervals and decide for the smartest split*/
SplitDecision LinearizeKKTWithMINRES::applySplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  std::vector<SplitApproximation> approximates(n_i_);
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
  abs_pos += n_st_x_;
  SmartPtr<const DenseVector> rhs_x_shift_upart = expandVector(rhs_shift_h_,rhs_dim_,abs_pos);
  abs_pos += n_sh_x_;
  SmartPtr<const DenseVector> rhs_u_part = expandVector(rhs_u_,rhs_dim_,abs_pos);
  abs_pos += n_u_;
  SmartPtr<const DenseVector> rhs_c_static_part = expandVector(rhs_static_c_,rhs_dim_,abs_pos);
  abs_pos += n_st_c_;
  SmartPtr<const DenseVector> rhs_c_shift_upart = expandVector(rhs_shift_c_,rhs_dim_,abs_pos);
  abs_pos += n_sh_c_;
  SmartPtr<const DenseVector> rhs_d_static_part = expandVector(rhs_static_d_,rhs_dim_,abs_pos);
  abs_pos += n_st_d_;
  SmartPtr<const DenseVector> rhs_d_shift_upart = expandVector(rhs_shift_d_,rhs_dim_,abs_pos);
  abs_pos += n_sh_d_;
  SmartPtr<const DenseVector> rhs_x_shift_lpart = expandVector(rhs_shift_h_,rhs_dim_,abs_pos);
  abs_pos += n_sh_x_;
  SmartPtr<const DenseVector> rhs_c_shift_lpart = expandVector(rhs_shift_c_,rhs_dim_,abs_pos);
  abs_pos += n_sh_c_;
  SmartPtr<const DenseVector> rhs_d_shift_lpart = expandVector(rhs_shift_d_,rhs_dim_,abs_pos);
  abs_pos += n_sh_d_;

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

/* expand Vector original from smaller (old) dim to large_dim, with the vector indices listing at which indices in the new vector the old values are to be found */  // disabled due to overloading issues
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

SmartPtr<const DenseVector> LinearizeKKTWithMINRES::expandVector(SmartPtr<const Vector> original,const std::vector<Index>& indices, const Index& large_dim) const
{
  const Index small_dim = indices.size();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d",small_dim,large_dim);
  assert(large_dim>=small_dim);

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

/*expand Vector original from smaller (old) dim to large_dim, inserting the original values as a dense block at start_idx in the new vector*/
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::expandVector(SmartPtr<const DenseVector> original,const Index& large_dim, const Index& start_idx) const
{
  return expandVector(dynamic_cast<const Vector*>(GetRawPtr(original)),large_dim,start_idx);
}


/*expand Vector original from smaller (old) dim to large_dim, inserting the original values as a dense block at start_idx in the new vector*/
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::expandVector(SmartPtr<const Vector> original,const Index& large_dim, const Index& start_idx) const
{
  const Index small_dim = original->Dim();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d   start_idx = %d",small_dim,large_dim,start_idx);
  assert(large_dim>=small_dim);

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


/* shrink Vector original from larger (old) dim to smaller dim, with the vector indices listing which indices from the old vector to keep*/
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::shrinkVector(SmartPtr<const DenseVector> original, const std::vector<Index>& indices) const
{
  const Index large_dim = original->Dim();
  const Index small_dim = indices.size();
//  printf("\nshrinkVector(): large_dim = %d   small_dim = %d",large_dim,small_dim);
  assert(large_dim>=small_dim);
  if (large_dim==small_dim) {
    //    printf("\nLinearizeKKTWithMINRES::shrinkVector(const): WARNING: called to shrink a vector to original size.");
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

/* shrink Vector original from larger (old) dim to smaller dim, with the vector indices listing which indices from the old vector to keep*/
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::shrinkVector(SmartPtr<const Vector> original, const std::vector<Index>& indices) const
{
  const Index large_dim = original->Dim();
  const Index small_dim = indices.size();
  //  printf("\nshrinkVector(): large_dim = %d   small_dim = %d",large_dim,small_dim);
  assert(large_dim>=small_dim);

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
  SmartPtr<DenseVector> new_vec = expanded_space->MakeNewDenseVector();

  em->TransMultVector(1.0,*original,0.0,*new_vec);

  return dynamic_cast<const DenseVector*>(GetRawPtr(new_vec));

}

/* same as shrinkVector, just with nonconst DenseVector entering*/
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::shrinkVector(SmartPtr<Vector> original, const std::vector<Index>& indices) const
{
  SmartPtr<const Vector> c_original = dynamic_cast<const Vector*>(GetRawPtr(original));

  return shrinkVector(c_original,indices);
}

Number LinearizeKKTWithMINRES::computeRMultS(SmartPtr<IpoptApplication> app) const
{
  Number retval;
  SmartPtr<const DenseVector> S = computeS(app);
  retval = rhs_i_->Dot(*S);

  return retval;
}

Number LinearizeKKTWithMINRES::computeSMultS(SmartPtr<IpoptApplication> app) const
{
  Number retval;
  Number squareval;

  SmartPtr<const DenseVector> S = computeS(app);
  retval = S->Dot(*S);

  return retval;
}

/* SmartPtr<const DenseVector> LinearizeKKTWithMINRES::computeR() const
{
// not needed for x_0 = 0 and 0 steps.
} */

SmartPtr<const DenseVector> LinearizeKKTWithMINRES::computeS(SmartPtr<IpoptApplication> app) const
{
  //initialize and set values where no further calculations are needed
  // static non-control parts of S equal static rhs part
  SmartPtr<const DenseVector> retval_st_h = expandVector(rhs_static_h_,rhs_dim_,0);
  SmartPtr<const DenseVector> retval_st_c = expandVector(rhs_static_c_,rhs_dim_,n_x_+n_u_);
  SmartPtr<const DenseVector> retval_st_d = expandVector(rhs_static_d_,rhs_dim_,n_x_+n_u_+n_c_);

  // init u part
  SmartPtr<DenseVector> ret_nc_u = dynamic_cast<DenseVector*>(rhs_u_->MakeNewCopy());
  SmartPtr<const DenseVector> retval_u = dynamic_cast<const DenseVector*>(rhs_u_->MakeNew());

  // upper shifted parts of S equal shift rhs part
  SmartPtr<const DenseVector> retval_ush_h = expandVector(rhs_shift_h_,rhs_dim_,n_st_x_);
  SmartPtr<const DenseVector> retval_ush_c = expandVector(rhs_shift_c_,rhs_dim_,n_x_+n_u_+n_st_c_);
  SmartPtr<const DenseVector> retval_ush_d = expandVector(rhs_shift_d_,rhs_dim_,n_x_+n_u_+n_c_+n_st_d_);

  // init lower shifted h part assigned to x
  SmartPtr<DenseVector> ret_nc_lsh_h = dynamic_cast<DenseVector*>(rhs_shift_h_->MakeNew());
  SmartPtr<const DenseVector> retval_lsh_h = dynamic_cast<const DenseVector*>(rhs_shift_h_->MakeNew());
  // init lower shifted c part
  SmartPtr<DenseVector> ret_nc_lsh_c = dynamic_cast<DenseVector*>(rhs_shift_c_->MakeNew());
  SmartPtr<const DenseVector> retval_lsh_c = dynamic_cast<const DenseVector*>(rhs_shift_c_->MakeNew());
  // init lower shifted d part
  SmartPtr<DenseVector> ret_nc_lsh_d = dynamic_cast<DenseVector*>(rhs_shift_d_->MakeNew());
  SmartPtr<const DenseVector> retval_lsh_d = dynamic_cast<const DenseVector*>(rhs_shift_d_->MakeNew());

  // get (lhs_i,shift times rhs_i) with i= (shift_h, u)
  SmartPtr<const DenseVector> tmp_rhs = expandVector(rhs_shift_h_,shift_x_indices_,lhs_h_->NCols());
  SmartPtr<DenseVector> extractor = dynamic_cast<DenseVector*>(tmp_rhs->MakeNew());
  lhs_h_->MultVector(1.0,*tmp_rhs,0.0,*extractor);

  // map result back to get u and shift part seperately for the appropriate parts of S
  SmartPtr<const DenseVector> extr_u_part = shrinkVector(extractor,u_indices_);
  ret_nc_u->AddOneVector(1.0,*extr_u_part,1.0);
  SmartPtr<const DenseVector> extr_shift_part = shrinkVector(extractor,shift_x_indices_);
  ret_nc_lsh_h->AddOneVector(1.0,*extr_shift_part,0.0);

  //get (lhs_c_shift^T times rhs_i) parts with i= (u, shift_c)
  tmp_rhs = expandVector(rhs_shift_c_,shift_c_indices_,lhs_c_->NRows());
  SmartPtr<DenseVectorSpace> extractor_space = new DenseVectorSpace(lhs_c_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  /*
  printf("\nlhs_c_->NCols() = %d, lhs_c_->NRows() = %d, tmp_rhs->Dim() = %d, extractor->Dim() = %d",lhs_c_->NCols(),lhs_c_->NRows(),tmp_rhs->Dim(),extractor->Dim());
  printf("\n\n");
  tmp_rhs->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "tmp_rhs");
  printf("\n\n");
  */
  lhs_c_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);

  //map results back to get u and shift part seperately for the appropriate parts of S
  extr_u_part = shrinkVector(extractor,u_indices_);
  ret_nc_u->AddOneVector(1.0,*extr_u_part,1.0);
  extr_shift_part = shrinkVector(extractor,shift_x_indices_);
  ret_nc_lsh_h->AddOneVector(1.0,*extr_shift_part,1.0);

  //get (lhs_d_shift^T times rhs_i) parts with i= (u, shift_d)
  tmp_rhs = expandVector(rhs_shift_d_,shift_d_indices_,lhs_d_->NRows());
  extractor_space = new DenseVectorSpace(lhs_d_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  lhs_d_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);

  //map result back to get u and shift part seperately for the appropriate parts of S
  extr_u_part = shrinkVector(extractor,u_indices_);
  ret_nc_u->AddOneVector(1.0,*extr_u_part,1.0);
  extr_shift_part = shrinkVector(extractor,shift_x_indices_);
  ret_nc_lsh_h->AddOneVector(1.0,*extr_shift_part,1.0);

  //map u part towards rhs dimension
  retval_u = dynamic_cast<const DenseVector*>(GetRawPtr(ret_nc_u));
  retval_u = expandVector(retval_u,rhs_dim_,n_x_);

  //get (lhs_c_shift times rhs_shift_c) part
  tmp_rhs = expandVector(rhs_shift_h_,shift_x_indices_,lhs_c_->NCols());
  extractor = dynamic_cast<DenseVector*>(tmp_rhs->MakeNew());
  lhs_c_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map result back to the appropriate parts of S
  extr_shift_part = shrinkVector(extractor,shift_c_indices_);
  ret_nc_lsh_c->AddOneVector(1.0,*extr_shift_part,0.0);

  //get (lhs_d_shift times rhs_shift_d) part
  tmp_rhs = expandVector(rhs_shift_h_,shift_x_indices_,lhs_d_->NCols());
  extractor = dynamic_cast<DenseVector*>(tmp_rhs->MakeNew());
  lhs_d_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map result back to the appropriate parts of S
  extr_shift_part = shrinkVector(extractor,shift_d_indices_);
  ret_nc_lsh_d->AddOneVector(1.0,*extr_shift_part,0.0);

  //map lower shifted parts towards rhs dimension
  retval_lsh_h = dynamic_cast<const DenseVector*>(GetRawPtr(ret_nc_lsh_h));
  retval_lsh_h = expandVector(retval_lsh_h,rhs_dim_,rhs_dim_-n_sh_x_-n_sh_c_-n_sh_d_);
  retval_lsh_c = dynamic_cast<const DenseVector*>(GetRawPtr(ret_nc_lsh_c));
  retval_lsh_c = expandVector(retval_lsh_c,rhs_dim_,rhs_dim_-n_sh_c_-n_sh_d_);
  retval_lsh_d = dynamic_cast<const DenseVector*>(GetRawPtr(ret_nc_lsh_d));
  retval_lsh_d = expandVector(retval_lsh_d,rhs_dim_,rhs_dim_-n_sh_d_);

  // get D_y_orig part
  extr_shift_part = computeDMultYrhs(app);

  // sum all the parts to get S
  SmartPtr<DenseVector> retval = extr_shift_part->MakeNewDenseVector();

  retval->AddTwoVectors(1.0,*retval_st_h,1.0,*retval_st_c,0.0);
  retval->AddTwoVectors(1.0,*retval_st_d,1.0,*retval_u,1.0);
  retval->AddTwoVectors(1.0,*retval_ush_h,1.0,*retval_ush_c,1.0);
  retval->AddTwoVectors(1.0,*retval_ush_d,1.0,*retval_lsh_h,1.0);
  retval->AddTwoVectors(1.0,*retval_lsh_c,1.0,*retval_lsh_d,1.0);
  retval->AddOneVector(1.0,*extr_shift_part,1.0);

  return dynamic_cast<const DenseVector*>(GetRawPtr(retval));
}
SmartPtr<const DenseVector> LinearizeKKTWithMINRES::computeDMultYrhs(SmartPtr<IpoptApplication> app) const
{
  // get KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

  // set up iterates vector and initialize - will be rhs for linear system
  SmartPtr<IteratesVector> it_vec = app->IpoptDataObject()->curr()->MakeNewIteratesVector();
  it_vec->Set(0.0);

  // set up static_h,shift_h and u parts of rhs iterates vector
  SmartPtr<DenseVectorSpace> orig_x_space = new DenseVectorSpace(n_x_+n_u_);
  SmartPtr<Vector> orig_x = orig_x_space->MakeNew();
  SmartPtr<const DenseVector> o_static = expandVector(rhs_static_h_,static_x_indices_,n_x_+n_u_);
  SmartPtr<const DenseVector> o_shift = expandVector(rhs_shift_h_,shift_x_indices_,n_x_+n_u_);
  // fill static and shift values into orig_x
  orig_x->AddTwoVectors(1.0,*o_static,1.0,*o_shift,0.0);
  // add control values
  o_static = expandVector(rhs_u_,u_indices_,n_x_+n_u_);
  orig_x->AddOneVector(1.0,*o_static,1.0);
  // set values into iterates vector
  it_vec->Set_x_NonConst(*orig_x);

  // set up static_c and shift_c parts of rhs iterates vector
  SmartPtr<DenseVectorSpace> orig_c_space = new DenseVectorSpace(n_c_);
  SmartPtr<Vector> orig_c = orig_c_space->MakeNew();
  o_static = expandVector(rhs_static_c_,static_c_indices_,n_c_);
  o_shift = expandVector(rhs_shift_c_,shift_c_indices_,n_c_);
  // fill static and shift values into orig_c
  orig_c->AddTwoVectors(1.0,*o_static,1.0,*o_shift,0.0);
  // set values into iterates vector
  it_vec->Set_y_c_NonConst(*orig_c);

  // set up static_d and shift_d parts of rhs iterates vector
  SmartPtr<DenseVectorSpace> orig_d_space = new DenseVectorSpace(n_d_);
  SmartPtr<Vector> orig_d = orig_d_space->MakeNew();
  o_static = expandVector(rhs_static_d_,static_d_indices_,n_d_);
  o_shift = expandVector(rhs_shift_d_,shift_d_indices_,n_d_);
  // fill static and shift values into orig_d
  orig_d->AddTwoVectors(1.0,*o_static,1.0,*o_shift,0.0);
  // set values into iterates vector
  it_vec->Set_y_d_NonConst(*orig_d);

  // do actual backsolve
  SmartPtr<IteratesVector> preretval = it_vec->MakeNewIteratesVector();
  pd_solver->Solve(1.0, 0.0, *it_vec, *preretval);

  // get (K^-1 times rhs) u part and eliminate non-u part elements
  SmartPtr<const Vector> retval_x = preretval->x();
  retval_x = dynamic_cast<const Vector*>(GetRawPtr(shrinkVector(retval_x,u_indices_)));
  retval_x = dynamic_cast<const Vector*>(GetRawPtr(expandVector(retval_x,u_indices_,lhs_h_->NCols())));
  SmartPtr<const Vector> retval_c;  // = preretval->y_c();
  SmartPtr<const Vector> retval_d;  // = preretval->y_d();

  // multiply with appropriate lhs parts
  SmartPtr<DenseVectorSpace> retval_nc_space = new DenseVectorSpace(lhs_h_->NRows());
  SmartPtr<Vector> retval_x_nonconst = retval_nc_space->MakeNew();
  lhs_h_->MultVector(1.0,*retval_x,0.0,*retval_x_nonconst);
  retval_nc_space = new DenseVectorSpace(lhs_c_->NRows());
  SmartPtr<Vector> retval_c_nonconst = retval_nc_space->MakeNew();
  lhs_c_->MultVector(1.0,*retval_x,0.0,*retval_c_nonconst);
  retval_nc_space = new DenseVectorSpace(lhs_d_->NRows());
  SmartPtr<Vector> retval_d_nonconst = retval_nc_space->MakeNew();
  lhs_d_->MultVector(1.0,*retval_x,0.0,*retval_d_nonconst);

  // reduce to the desired quantities
  retval_x = dynamic_cast<const Vector*>(GetRawPtr(shrinkVector(retval_x_nonconst,shift_x_indices_)));
  retval_c = dynamic_cast<const Vector*>(GetRawPtr(shrinkVector(retval_c_nonconst,shift_c_indices_)));
  retval_d = dynamic_cast<const Vector*>(GetRawPtr(shrinkVector(retval_d_nonconst,shift_d_indices_)));

  //map upwards to rhs dimension
  Index abs_pos = rhs_dim_-retval_x->Dim()-retval_c->Dim()-retval_d->Dim();
  retval_x = dynamic_cast<const Vector*>(GetRawPtr(expandVector(retval_x,rhs_dim_,abs_pos)));
  abs_pos += n_sh_x_;
  retval_c = dynamic_cast<const Vector*>(GetRawPtr(expandVector(retval_c,rhs_dim_,abs_pos)));
  abs_pos += n_sh_c_;
  retval_d = dynamic_cast<const Vector*>(GetRawPtr(expandVector(retval_d,rhs_dim_,abs_pos)));

  // add all results to get resulting D-part of S
  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(rhs_dim_);
  SmartPtr<DenseVector> retval = retval_space->MakeNewDenseVector();
  retval->AddTwoVectors(1.0,*retval_x,1.0,*retval_c,0.0);
  retval->AddOneVector(1.0,*retval_d,1.0);

  return dynamic_cast<const DenseVector*>(GetRawPtr(retval));
}

SmartPtr<const DenseVector> LinearizeKKTWithMINRES::undoConditionning(SmartPtr<IpoptApplication> app, SmartPtr<Vector> cond_sol) const
{
  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(rhs_dim_);

  // split static, control and upper shift parts from the conditionned solution into h,c and d segments
  std::vector<Index> st_h_indices(n_st_x_);
  for (int i=0;i<st_h_indices.size();i++)
    st_h_indices[i] = i;
  SmartPtr<const DenseVector> retval_st_h_part = shrinkVector(cond_sol,st_h_indices);

  std::vector<Index> ush_h_indices(n_sh_x_);
  for (int i=0;i<ush_h_indices.size();i++)
    ush_h_indices[i] = i+n_st_x_;
  SmartPtr<const DenseVector> retval_ush_h_part = shrinkVector(cond_sol,ush_h_indices);

  std::vector<Index> u_indices(n_u_);
  for (int i=0;i<u_indices.size();i++)
    u_indices[i] = i+n_x_;
  SmartPtr<const DenseVector> retval_u_part = shrinkVector(cond_sol,u_indices);

  std::vector<Index> st_c_indices(n_st_c_);
  for (int i=0;i<st_c_indices.size();i++)
    st_c_indices[i] = i+n_x_+n_u_;
  SmartPtr<const DenseVector> retval_st_c_part = shrinkVector(cond_sol,st_c_indices);

  std::vector<Index> ush_c_indices(n_sh_c_);
  for (int i=0;i<ush_c_indices.size();i++)
    ush_c_indices[i] = i+n_x_+n_u_+n_st_c_;
  SmartPtr<const DenseVector> retval_ush_c_part = shrinkVector(cond_sol,ush_c_indices);

  std::vector<Index> st_d_indices(n_st_d_);
  for (int i=0;i<st_d_indices.size();i++)
    st_d_indices[i] = i+n_x_+n_u_+n_c_;
  SmartPtr<const DenseVector> retval_st_d_part = shrinkVector(cond_sol,st_d_indices);

  std::vector<Index> ush_d_indices(n_sh_d_);
  for (int i=0;i<ush_d_indices.size();i++)
    ush_d_indices[i] = i+n_x_+n_u_+n_c_+n_st_d_;
  SmartPtr<const DenseVector> retval_ush_d_part = shrinkVector(cond_sol,ush_d_indices);

  // remap the parts to fit the original problem formulation
  retval_st_h_part = expandVector(retval_st_h_part,static_x_indices_,n_x_+n_u_);
  retval_ush_h_part = expandVector(retval_ush_h_part,shift_x_indices_,n_x_+n_u_);
  retval_u_part = expandVector(retval_u_part,u_indices_,n_x_+n_u_);
  retval_st_c_part = expandVector(retval_st_c_part,static_c_indices_,n_c_);
  retval_ush_c_part = expandVector(retval_ush_c_part,shift_c_indices_,n_c_);
  retval_st_d_part = expandVector(retval_st_d_part,static_d_indices_,n_d_);
  retval_ush_d_part = expandVector(retval_ush_d_part,shift_d_indices_,n_d_);

  // get KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

  // set up iterates vector and initialize - will be rhs for linear system
  SmartPtr<IteratesVector> it_vec = app->IpoptDataObject()->curr()->MakeNewIteratesVector();
  it_vec->Set(0.0);

  // set up static_h,shift_h and u parts of rhs iterates vector
  SmartPtr<Vector> orig_x = retval_st_h_part->MakeNew();
  // fill static and shift values into orig_x
  orig_x->AddTwoVectors(1.0,*retval_st_h_part,1.0,*retval_ush_h_part,0.0);
  // add control values
  orig_x->AddOneVector(1.0,*retval_u_part,1.0);
  // set values into iterates vector
  it_vec->Set_x_NonConst(*orig_x);

  // set up static_c and shift_c parts of rhs iterates vector
  SmartPtr<Vector> orig_c = retval_st_c_part->MakeNew();
  // fill static and shift values into orig_c
  orig_c->AddTwoVectors(1.0,*retval_st_c_part,1.0,*retval_ush_c_part,0.0);
  // set values into iterates vector
  it_vec->Set_y_c_NonConst(*orig_c);

  // set up static_d and shift_d parts of rhs iterates vector
  SmartPtr<Vector> orig_d = retval_st_d_part->MakeNew();
  // fill static and shift values into orig_d
  orig_d->AddTwoVectors(1.0,*retval_st_d_part,1.0,*retval_ush_d_part,0.0);
  // set values into iterates vector
  it_vec->Set_y_d_NonConst(*orig_d);

  // do actual backsolve
  SmartPtr<IteratesVector> preretval = it_vec->MakeNewIteratesVector();
  pd_solver->Solve(1.0, 0.0, *it_vec, *preretval);

  SmartPtr<const Vector> retval_x = preretval->x();
  SmartPtr<const Vector> retval_c = preretval->y_c();
  SmartPtr<const Vector> retval_d = preretval->y_d();

  // extract results
  retval_st_h_part = shrinkVector(retval_x,static_x_indices_);
  retval_ush_h_part = shrinkVector(retval_x,shift_x_indices_);
  retval_u_part = shrinkVector(retval_x,u_indices_);
  retval_st_c_part = shrinkVector(retval_c,static_c_indices_);
  retval_ush_c_part = shrinkVector(retval_c,shift_c_indices_);
  retval_st_d_part = shrinkVector(retval_d,static_d_indices_);
  retval_ush_d_part = shrinkVector(retval_d,shift_d_indices_);

  // expand to final dimension
  Index abs_pos = 0;
  retval_st_h_part = expandVector(retval_st_h_part,rhs_dim_,abs_pos);
  abs_pos += n_st_x_;
  retval_ush_h_part = expandVector(retval_ush_h_part,rhs_dim_,abs_pos);
  abs_pos += n_sh_x_;
  retval_u_part = expandVector(retval_u_part,rhs_dim_,abs_pos);
  abs_pos += n_u_;
  retval_st_c_part = expandVector(retval_st_c_part,rhs_dim_,abs_pos);
  abs_pos += n_st_c_;
  retval_ush_c_part = expandVector(rhs_shift_c_,rhs_dim_,abs_pos);
  abs_pos += n_sh_c_;
  retval_st_d_part = expandVector(retval_st_d_part,rhs_dim_,abs_pos);
  abs_pos += n_st_d_;
  retval_ush_d_part = expandVector(retval_ush_d_part,rhs_dim_,abs_pos);

  // copy lower shifted part over to retval part
  std::vector<Index> lsh_indices(n_sh_x_+n_sh_c_+n_sh_d_);
  for (int i=0;i<lsh_indices.size();i++)
    lsh_indices[i] = n_x_+n_u_+n_c_+n_d_+i;
  SmartPtr<const DenseVector> retval_lsh_part = shrinkVector(cond_sol,lsh_indices);
  retval_lsh_part = expandVector(dynamic_cast<const Vector*>(GetRawPtr(retval_lsh_part)),lsh_indices,rhs_dim_);

  SmartPtr<DenseVector> retval = retval_space->MakeNewDenseVector();
      /*      printf("\n");
      for (int i=0;i<n_x_+n_c_+n_d_+n_u_;i++)
	printf("\n different");

retval_st_h_part->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"retval_st_h_part");
      printf("\n\n");
retval_ush_h_part->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"retval_ush_h_part");
      printf("\n\n");
retval_u_part->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"retval_u_part");
      printf("\n\n");
retval_st_c_part->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"retval_st_c_part");
      printf("\n\n");
retval_ush_c_part->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"retval_ush_c_part");
      printf("\n\n");
retval_st_d_part->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"retval_st_d_part");
      printf("\n\n");
retval_ush_d_part->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"retval_ush_d_part");
            printf("\n\n");
retval_lsh_part->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"retval_lsh_part");
      printf("\n\n");
*/
  retval->AddTwoVectors(1.0,*retval_st_h_part,1.0,*retval_ush_h_part,0.0);
  retval->AddTwoVectors(1.0,*retval_u_part,1.0,*retval_st_c_part,1.0);
  retval->AddTwoVectors(1.0,*retval_ush_c_part,1.0,*retval_st_d_part,1.0);
  retval->AddTwoVectors(1.0,*retval_ush_d_part,1.0,*retval_lsh_part,1.0);

  return dynamic_cast<const DenseVector*>(GetRawPtr(retval));
}

SmartPtr<const DenseVector> LinearizeKKTWithMINRES::computeDeltaP(const IntervalInfoSet& intervals,const Index& interval) const
{
  const Index n_parameters = intervals.getParameterCount();
  Index counter = 0;
  std::vector<Index> parameterflags = intervals.getParameterIDVec();
  std::vector<Index> intervalflags = intervals.getIntervalIDVec();
  std::vector<Index> parameterIDs = intervals.getParameterIDs();
  IntWidthScaling *int_widthscaler = new IntWidthScaling();
  std::vector<Number> int_widths = int_widthscaler->scaleIntervalWidths(intervals,interval);
  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(2*n_parameters);
  SmartPtr<DenseVector> retval = retval_space->MakeNewDenseVector();
  Number* ret_values = new Number[2*n_parameters];

  for (int k=0; k<intervalflags.size(); k++) {
    for (int j=0; j<parameterIDs.size(); j++) {
      if (intervalflags[k] == interval) {
	if (parameterflags[k] == parameterIDs[j]) {
	  if (intervals.isUpper(k))
	    ret_values[counter] = -0.5*int_widths[j];
	   else
	    ret_values[counter] = 0.5*int_widths[j];
	  counter++;
	}
      }
    }
  }

  retval->SetValues(ret_values);
  return dynamic_cast<const DenseVector*>(GetRawPtr(retval));
}

SplitDecision LinearizeKKTWithMINRES::decideInterval(const std::vector<SplitApproximation>& approximates) const
{
  /* waiting to be implemeted*/
}

SplitApproximation LinearizeKKTWithMINRES::applyAlgorithmOnInterval(SmartPtr<IpoptApplication> app, const Index& interval)
{
  // assign private attributes depending on chosen interval
  std::vector<Index> intervalIDs = intervals_.getIntervalIDVec();
  const Index n_parameters = intervals_.getParameterCount();
  SmartPtr<DenseVectorSpace> z_space = new DenseVectorSpace(n_x_+n_u_+n_c_+n_sh_x_+n_sh_c_+n_sh_d_);
  SmartPtr<Vector> z_postsplit_cond = z_space->MakeNewDenseVector();
  SmartPtr<const DenseVector> z_postsplit = z_space->MakeNewDenseVector();
  SmartPtr<const DenseVector> deltap;
  SmartPtr<MultiVectorMatrixSpace> sense_space = new MultiVectorMatrixSpace(n_parameters*2,*z_space);
  SmartPtr<MultiVectorMatrix> sense = sense_space->MakeNewMultiVectorMatrix();
  SmartPtr<DenseVectorSpace> res_space = new DenseVectorSpace(sense->NRows());
  SmartPtr<DenseVector> res_sense = res_space->MakeNewDenseVector();
  Index s_count = 0;
  for (Index j=0;j<n_p_;j++) {
    if (intervalIDs[j]==interval) {
      if (assignIntAndParaDepAttr(interval,j)) {
	// compute the solution of the MINRES-Step
	const Number steplength = computeRMultS(app)*computeSMultS(app);
	z_postsplit_cond = rhs_i_->MakeNewCopy();
	// printf("\n");
	// z_postsplit_cond->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"z_slow");
	// printf("\n");
	if (steplength) {
	  //	printf("\nLinearizeKKTWithMINRES::applyAlgorithmOnInterval(): MINRES-Step is %e.",steplength);
	  z_postsplit_cond->Scal(steplength);
	  //	  z_postsplit_cond->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"z_sc_slow");
	  z_postsplit = computeS(app);
	  z_postsplit->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"S_slow");
	} else
	  printf("\nLinearizeKKTWithMINRES::applyAlgorithmOnInterval(): ERROR: MINRES-Step is 0.");

	// undo conditionning step to get the approximation of the postsplit sensitivity column
	z_postsplit = undoConditionning(app,z_postsplit_cond);

	// compare sense vectors for test purposes
	//      SmartPtr<Vector> compy = z_postsplit->MakeNew();
	/*      printf("\n");
		for (int i=0;i<n_x_+n_c_+n_d_+n_u_;i++)
		printf("\n different");
		//      compy->AddTwoVectors(1.0,*z_copy,-1.0,*z_postsplit_cond,0.0);
		printf("\n\n");
		z_postsplit->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"z_post");
		//      compy->AddTwoVectors(1.0,*z_copy,-1.0,*z_postsplit,0.0);
		printf("\n\n");
		//      compy->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"compycomp");

		*/

	// add the resulting sense vector to matrix ....
	sense->SetVector(s_count,*z_postsplit);
	s_count++;


      } else
	printf("\nLinearizeKKTWithMINRES::applyAlgorithmOnInterval(): ERROR: Unable to assign interval dependent attributes!");
    }
  }
  // get the parameter perturbation for this split
  deltap = computeDeltaP(intervals_,interval);

  // calculate absolute sensitivity values
  sense->MultVector(1.0,*deltap,0.0,*res_sense);

  printf("\n\n");
  res_sense->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"result_s");
  printf("\n\n");
  std::vector<Index> x_intervalIDs = (dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_->OwnerSpace())))->GetIntegerMetaData("intervalID");
  printf("\n");

  // for (int i=0;i<x_intervalIDs.size();i++)
  //       printf("\n interval---%d----",x_intervalIDs[i]);
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
  SmartPtr<OptionsList> options = app->Options();
  // read from options what kind of controlselector to return
  std::string sensemode;
  SplitAlgorithm* retval = new LinearizeKKTWithMINRES(app);
  if (options->GetStringValue("sensemode",sensemode,"")) {
    if (sensemode == "MINRES") {
      retval =  new LinearizeKKTWithMINRES(app);
    } if (sensemode == "control") {
      retval = new SplitWRTControlSensitivities();
    } if (sensemode == "GMRES") {
      retval = new LinKKTFaster(app);
    }
  } else
    printf("\nassignSplitAlgorithm(): ERROR: No splitAlgorithm chosen!");
  return retval;
}

LinKKTFaster::LinKKTFaster(SmartPtr<IpoptApplication> app)
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
  intervals_ = IntervalInfoSet(p_);
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
    //   printf("LinKKTFaster::LinKKTFaster(): c_intervalIDs_[%d]=%d\n",i,c_intervalIDs_[i]);
    // }
  } else
    printf("LinKKTFaster::LinKKTFaster(): Error - no intervalIDs for c_space!\n");
  if (d_space->HasIntegerMetaData("intervalID")) {
    d_intervalIDs_ = d_space->GetIntegerMetaData("intervalID");
    // for (int i=0;i<d_intervalIDs_.size();i++) {
    //   printf("LinKKTFaster::LinKKTFaster(): d_intervalIDs_[%d]=%d\n",i,d_intervalIDs_[i]);
    // }
  } else
    printf("LinKKTFaster::LinKKTFaster(): Error - no intervalIDs for d_space!\n");

  n_i_ = intervals_.getIntervalCount();
  n_u_ = u_indices_.size();
  n_p_ = intervals_.Size();
  n_x_ = x_intervalIDs_.size()-n_u_;
  n_c_ = c_intervalIDs_.size();
  n_d_ = d_intervalIDs_.size();

  n_sh_x_ = int(n_x_/n_i_);
  n_st_x_ = n_x_ - n_sh_x_;
  n_sh_c_ = int(n_c_/n_i_);
  n_st_c_ = n_c_ - n_sh_c_;
  n_sh_d_ = int(n_d_/n_i_);
  n_st_d_ = n_d_ - n_sh_d_;

  rhs_dim_ = computeRHSDim();

  // interval specific extracted data cannot be initialized in constructor!!

}

/* apply MINRES-approximation of a split to all intervals and decide for the smartest split*/
SplitDecision LinKKTFaster::applySplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  std::vector<SplitApproximation> approximates(n_i_);
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

Index LinKKTFaster::computeRHSDim() const
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
  //  printf("Index LinKKTFaster::computeRHSDim(): retval = %d,   n_x_ = %d, n_c_ = %d, n_d_ = %d,   n_u_ = %d ",retval, n_x_, n_c_, n_d_, n_u_);
  retval += n_x_ +n_c_ +n_d_+n_u_;

  return retval;
}

bool LinKKTFaster::assignIntAndParaDepAttr(const Index& interval,const Index& column)
{
  bool retval = 1;

    // split interval indices into to be shifted interval entries and remainder
  if (!splitIntervalIndices(x_intervalIDs_,static_x_indices_,shift_x_indices_,interval)) {
    printf("\nLinKKTFaster::assignIntAndParaDepAttr(): ERROR: unable to split x_intervalIDs_");
    retval = 0;
  }
  /*  for (int i=0;i<static_x_indices_.size();i++)
      printf("\nstatic_x_indices_[%d] = %d",i,static_x_indices_[i]);
  printf("\n\n");

  */

  if (!splitIntervalIndices(c_intervalIDs_,static_c_indices_,shift_c_indices_,interval)) {
    printf("\nLinKKTFaster::assignIntAndParaDepAttr(): ERROR: unable to split c_intervalIDs_");
    retval = 0;
  }
  if (!splitIntervalIndices(d_intervalIDs_,static_d_indices_,shift_d_indices_,interval)) {
    printf("\nLinKKTFaster::assignIntAndParaDepAttr(): ERROR: unable to split d_intervalIDs_");
    retval = 0;
  }
  // get different rhs entries for this column
  SmartPtr<const DenseVector> rhs_h_i;
  if (GetRawPtr(rhs_h_)) {
    x_i_ = extractColumn(rhs_h_,column);
    rhs_h_i = x_i_;
  }
  else
    printf("\nLinKKTFaster::assignIntAndParaDepAttr(): ERROR: rhs_h_ is NULL");
  SmartPtr<const DenseVector> rhs_c_i;
  if (GetRawPtr(rhs_c_)) {
    y_c_i_ = extractColumn(rhs_c_,column);
    rhs_c_i = y_c_i_;
  }
  else
    printf("\nLinKKTFaster::assignIntAndParaDepAttr(): ERROR: rhs_c_ is NULL");
  SmartPtr<const DenseVector> rhs_d_i;
  if (GetRawPtr(rhs_d_)) {
    rhs_d_i = extractColumn(rhs_d_,column);
    y_d_i_ = rhs_d_i;
  }
  else
    printf("\nLinKKTFaster::assignIntAndParaDepAttr(): ERROR: rhs_d_ is NULL");



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
  abs_pos += n_st_x_;
  SmartPtr<const DenseVector> rhs_x_shift_upart = expandVector(rhs_shift_h_,rhs_dim_,abs_pos);
  abs_pos += n_sh_x_;
  SmartPtr<const DenseVector> rhs_u_part = expandVector(rhs_u_,rhs_dim_,abs_pos);
  abs_pos += n_u_;
  SmartPtr<const DenseVector> rhs_c_static_part = expandVector(rhs_static_c_,rhs_dim_,abs_pos);
  abs_pos += n_st_c_;
  SmartPtr<const DenseVector> rhs_c_shift_upart = expandVector(rhs_shift_c_,rhs_dim_,abs_pos);
  abs_pos += n_sh_c_;
  SmartPtr<const DenseVector> rhs_d_static_part = expandVector(rhs_static_d_,rhs_dim_,abs_pos);
  abs_pos += n_st_d_;
  SmartPtr<const DenseVector> rhs_d_shift_upart = expandVector(rhs_shift_d_,rhs_dim_,abs_pos);
  abs_pos += n_sh_d_;
  SmartPtr<const DenseVector> rhs_x_shift_lpart = expandVector(rhs_shift_h_,rhs_dim_,abs_pos);
  abs_pos += n_sh_x_;
  SmartPtr<const DenseVector> rhs_c_shift_lpart = expandVector(rhs_shift_c_,rhs_dim_,abs_pos);
  abs_pos += n_sh_c_;
  SmartPtr<const DenseVector> rhs_d_shift_lpart = expandVector(rhs_shift_d_,rhs_dim_,abs_pos);
  abs_pos += n_sh_d_;

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
SmartPtr<const DenseVector> LinKKTFaster::extractColumn(SmartPtr<const Matrix> original,const Index& column) const
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

/* expand Vector original from smaller (old) dim to large_dim, with the vector indices listing at which indices in the new vector the old values are to be found */  // disabled due to overloading issues
SmartPtr<const DenseVector> LinKKTFaster::expandVector(SmartPtr<const DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const
{
  const Index small_dim = indices.size();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d",small_dim,large_dim);
  assert(large_dim>=small_dim);

  if (large_dim==small_dim) {
//    printf("\nLinKKTFaster::expandVector(): WARNING: called to expand a vector to original size.");
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
  printf("\nLinKKTFaster::expandVector(): Unknown ERROR.");
  return NULL;
}

SmartPtr<const DenseVector> LinKKTFaster::expandVector(SmartPtr<DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const
{
  SmartPtr<const DenseVector> c_original = dynamic_cast<const DenseVector*>(GetRawPtr(original));
  return expandVector(c_original,indices,large_dim);
}

SmartPtr<const DenseVector> LinKKTFaster::expandVector(SmartPtr<const Vector> original,const std::vector<Index>& indices, const Index& large_dim) const
{
  const Index small_dim = indices.size();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d",small_dim,large_dim);
  assert(large_dim>=small_dim);

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

/*expand Vector original from smaller (old) dim to large_dim, inserting the original values as a dense block at start_idx in the new vector*/
SmartPtr<const DenseVector> LinKKTFaster::expandVector(SmartPtr<const DenseVector> original,const Index& large_dim, const Index& start_idx) const
{
  return expandVector(dynamic_cast<const Vector*>(GetRawPtr(original)),large_dim,start_idx);
}


/*expand Vector original from smaller (old) dim to large_dim, inserting the original values as a dense block at start_idx in the new vector*/
SmartPtr<const DenseVector> LinKKTFaster::expandVector(SmartPtr<const Vector> original,const Index& large_dim, const Index& start_idx) const
{
  const Index small_dim = original->Dim();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d   start_idx = %d",small_dim,large_dim,start_idx);
  assert(large_dim>=small_dim);

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

SmartPtr<DenseVector> LinKKTFaster::expand(SmartPtr<DenseVector> original,const Index& large_dim, const Index& start_idx) const
{
  const Index small_dim = original->Dim();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d   start_idx = %d",small_dim,large_dim,start_idx);
  assert(large_dim>=small_dim);

  Index* exppos = new Index[small_dim];
  for (int i=0;i<small_dim;i++)
    exppos[i]=i+start_idx;
  //    printf("\n");
  //    for (int i=0;i<small_dim;i++)
  //      printf("\nexppos[%d] = %d",i,exppos[i]);

  SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
  SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> expanded_space = new DenseVectorSpace(large_dim);
  SmartPtr<DenseVector> new_vec = expanded_space->MakeNewDenseVector();

  em->MultVector(1.0,*original,0.0,*new_vec);

  return new_vec;
}

SmartPtr<DenseVector> LinKKTFaster::expand(SmartPtr<IteratesVector> original,const Index& large_dim, const Index& start_idx) const
{
  // extract relevant parts of IteratesVector and transform into DenseVectors
  const Index total_dim = original->x()->Dim()+original->y_c()->Dim()+original->y_d()->Dim();
  Index abs_pos = 0;
  SmartPtr<DenseVector> x_part = expand(dynamic_cast<DenseVector*>(original->x()->MakeNewCopy()),total_dim,abs_pos);
  abs_pos += original->x()->Dim();
  SmartPtr<DenseVector> y_c_part = expand(dynamic_cast<DenseVector*>(original->y_c()->MakeNewCopy()),total_dim,abs_pos);
  abs_pos += original->y_c()->Dim();
  SmartPtr<DenseVector> y_d_part = expand(dynamic_cast<DenseVector*>(original->y_d()->MakeNewCopy()),total_dim,abs_pos);
  abs_pos += original->y_d()->Dim();
  SmartPtr<DenseVectorSpace> ori_space = new DenseVectorSpace(total_dim);
  SmartPtr<DenseVector> new_ori = ori_space->MakeNewDenseVector();

  // sum DenseVectors
  new_ori->AddTwoVectors(1.0,*x_part,1.0,*y_c_part,0.0);
  new_ori->AddOneVector(1.0,*y_d_part,1.0);

  return expand(new_ori,large_dim,start_idx);
}

/*expand Vector original from smaller (old) dim to large_dim, inserting the original values as a dense block at start_idx in the new vector*/
SmartPtr<DenseVector> LinKKTFaster::expand(SmartPtr<DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const
{
  const Index small_dim = indices.size();
//  printf("\nexpandVector(): small_dim = %d   large_dim = %d",small_dim,large_dim);
  assert(large_dim>=small_dim);

  if (large_dim==small_dim) {
//    printf("\nLinKKTFaster::expandVector(): WARNING: called to expand a vector to original size.");
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
    SmartPtr<DenseVector> new_vec = expanded_space->MakeNewDenseVector();

    em->MultVector(1.0,*original,0.0,*new_vec);

    return new_vec;
  }
  printf("\nLinKKTFaster::expand(): Unknown ERROR.");
  return NULL;
}


/* shrink Vector original from larger (old) dim to smaller dim, with the vector indices listing which indices from the old vector to keep*/
SmartPtr<DenseVector> LinKKTFaster::shrink(SmartPtr<DenseVector> original,const std::vector<Index>& indices) const
{
  const Index large_dim = original->Dim();
  const Index small_dim = indices.size();
//  printf("\nshrinkVector(): large_dim = %d   small_dim = %d",large_dim,small_dim);
  assert(large_dim>=small_dim);
  if (large_dim==small_dim) {
    //    printf("\nLinKKTFaster::shrink(): WARNING: called to shrink a vector to original size.");
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
    SmartPtr<DenseVector> new_vec = expanded_space->MakeNewDenseVector();

    em->TransMultVector(1.0,*original,0.0,*new_vec);

    return new_vec;
  }
  printf("\nLinKKTFaster::shrink(): Unknown ERROR.");
  return NULL;
}

/* shrink Vector original from larger (old) dim to smaller dim, with the vector indices listing which indices from the old vector to keep*/
SmartPtr<const DenseVector> LinKKTFaster::shrinkVector(SmartPtr<const DenseVector> original, const std::vector<Index>& indices) const
{
  const Index large_dim = original->Dim();
  const Index small_dim = indices.size();
//  printf("\nshrinkVector(): large_dim = %d   small_dim = %d",large_dim,small_dim);
  assert(large_dim>=small_dim);
  if (large_dim==small_dim) {
    //    printf("\nLinKKTFaster::shrinkVector(const): WARNING: called to shrink a vector to original size.");
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
  printf("\nLinKKTFaster::shrinkVector(const): Unknown ERROR.");
  return NULL;
}

/* same as shrinkVector, just with nonconst DenseVector entering */
   SmartPtr<const DenseVector> LinKKTFaster::shrinkVector(SmartPtr<DenseVector> original, const std::vector<Index>& indices) const
{
  SmartPtr<const DenseVector> c_original = dynamic_cast<const DenseVector*>(GetRawPtr(original));

  return shrinkVector(c_original,indices);
}

/* shrink Vector original from larger (old) dim to smaller dim, with the vector indices listing which indices from the old vector to keep*/
SmartPtr<const DenseVector> LinKKTFaster::shrinkVector(SmartPtr<const Vector> original, const std::vector<Index>& indices) const
{
  const Index large_dim = original->Dim();
  const Index small_dim = indices.size();
  //  printf("\nshrinkVector(): large_dim = %d   small_dim = %d",large_dim,small_dim);
  assert(large_dim>=small_dim);

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
  SmartPtr<DenseVector> new_vec = expanded_space->MakeNewDenseVector();

  em->TransMultVector(1.0,*original,0.0,*new_vec);

  return dynamic_cast<const DenseVector*>(GetRawPtr(new_vec));

}

/* same as shrinkVector, just with nonconst DenseVector entering */
SmartPtr<const DenseVector> LinKKTFaster::shrinkVector(SmartPtr<Vector> original, const std::vector<Index>& indices) const
{
  SmartPtr<const Vector> c_original = dynamic_cast<const Vector*>(GetRawPtr(original));

  return shrinkVector(c_original,indices);
}

Number LinKKTFaster::computeRMultS(SmartPtr<IpoptApplication> app,const VectorSet& target) const
{
  Number retval =0;

  VectorSet S = computeS(app,target);
  retval+= target.top->Dot(*S.top)+target.x->Dot(*S.x)+target.y_c->Dot(*S.y_c)+target.y_d->Dot(*S.y_d);
  return retval;
}

Number LinKKTFaster::computeSMultS(SmartPtr<IpoptApplication> app,const VectorSet& target) const
{
  Number retval =0;

  VectorSet S = computeS(app,target);
  retval+= S.top->Dot(*S.top)+S.x->Dot(*S.x)+S.y_c->Dot(*S.y_c)+S.y_d->Dot(*S.y_d);

  return retval;
}

/* SmartPtr<const DenseVector> LinKKTFaster::computeR() const
{
// not needed for x_0 = 0 and 0 steps.
} */

VectorSet LinKKTFaster::computeS(SmartPtr<IpoptApplication> app,const VectorSet& target) const
{
  // init retval parts
  SmartPtr<DenseVector> retval_u = rhs_u_->MakeNewDenseVector();
  SmartPtr<DenseVector> retval_x = rhs_shift_h_->MakeNewDenseVector();
  SmartPtr<DenseVector> retval_y_c = rhs_shift_c_->MakeNewDenseVector();
  SmartPtr<DenseVector> retval_y_d = rhs_shift_d_->MakeNewDenseVector();

  // get (W_i,shift times rhs_i) with i= (shift_h, u)
  SmartPtr<DenseVector> tmp_rhs = expand(target.x,shift_x_indices_,lhs_h_->NCols());
  SmartPtr<DenseVector> extractor = tmp_rhs->MakeNewDenseVector();
  lhs_h_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  // map results back to get u and shift part seperately for the appropriate parts of S
  SmartPtr<DenseVector> extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,0.0);
  SmartPtr<DenseVector> extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,0.0);

  //get (A_shift^T times rhs_i) parts with i= (u, shift_c)
  tmp_rhs = expand(target.y_c,shift_c_indices_,lhs_c_->NRows());
  SmartPtr<DenseVectorSpace> extractor_space = new DenseVectorSpace(lhs_c_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  lhs_c_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map results back to get u and shift part seperately for the appropriate parts of S
  extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,1.0);
  extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,1.0);

  //get (B_shift^T times rhs_i) parts with i= (u, shift_d)
  tmp_rhs = expand(target.y_d,shift_d_indices_,lhs_d_->NRows());
  extractor_space = new DenseVectorSpace(lhs_d_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  lhs_d_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map result back to get u and shift part seperately for the appropriate parts of S
  extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,1.0);
  extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,1.0);

  //get (A_shift times rhs_shift_c) part
  tmp_rhs = expand(target.x,shift_x_indices_,lhs_c_->NCols());
  extractor = tmp_rhs->MakeNewDenseVector();
  lhs_c_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map result back to the appropriate parts of S
  retval_y_c = shrink(extractor,shift_c_indices_);

  //get (B_shift times rhs_shift_d) part
  tmp_rhs = expand(target.x,shift_x_indices_,lhs_d_->NCols());
  extractor = tmp_rhs->MakeNewDenseVector();
  lhs_d_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map result back to the appropriate parts of S
  retval_y_d = shrink(extractor,shift_d_indices_);

  // copy top values over from target, at the same time add D-multiplication to shifted part
  VectorSet retval = computeDMultYrhs(app,target);

  // map and insert u part back into top part
  SmartPtr<DenseVector> top_x = dynamic_cast<DenseVector*>(target.top->x()->MakeNewCopy());
  SmartPtr<const DenseVector> top_u = expandVector(retval_u,u_indices_,top_x->Dim());
  top_x->AddOneVector(1.0,*top_u,1.0);
  SmartPtr<IteratesVector> retval_top = target.top->MakeNewIteratesVector();
  retval_top->Set(0.0);
  retval_top->Set_x_NonConst(*top_x);
  retval_top->Set_y_c(*target.top->y_c());
  retval_top->Set_y_d(*target.top->y_d());

  // insert calculated parts into retval
  retval.top = retval_top;
  retval.x->AddOneVector(1.0,*retval_x,1.0);
  retval.y_c->AddOneVector(1.0,*retval_y_c,1.0);
  retval.y_d->AddOneVector(1.0,*retval_y_d,1.0);

  return retval;
}

  SmartPtr<DenseVector> LinKKTFaster::resetValuesAt(SmartPtr<DenseVector> target, std::vector<Index>indices) const
{
  assert(indices.size()<=target->Dim());
  const Index dim = target->Dim();
  const Number* values = target->Values();
  Number* new_vals = new Number[dim];
  std::copy(values,values+dim,new_vals);

  for (int i=0;i<indices.size();i++) {
    assert(indices[i]<target->Dim());
    new_vals[u_indices_[i]] = 0.0;
  }
  SmartPtr<DenseVector> retval = target->MakeNewDenseVector();
  retval->SetValues(new_vals);

  return retval;
}

VectorSet LinKKTFaster::computeDMultYrhs(SmartPtr<IpoptApplication> app,const VectorSet& target) const
{
  // get KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

  // do actual backsolve
  SmartPtr<IteratesVector> preretval = target.top->MakeNewIteratesVector();
  pd_solver->Solve(1.0, 0.0, *target.top, *preretval);

  // get (K^-1 times rhs) u part and eliminate non-u part elements
  SmartPtr<DenseVector> retval_x = dynamic_cast<DenseVector*>(preretval->x()->MakeNewCopy());
  retval_x = shrink(retval_x,u_indices_);
  retval_x = expand(retval_x,u_indices_,lhs_h_->NCols());
  SmartPtr<DenseVector> retval_c;  // = preretval->y_c();
  SmartPtr<DenseVector> retval_d;  // = preretval->y_d();

  // multiply with appropriate lhs parts
  SmartPtr<DenseVectorSpace> retval_nc_space = new DenseVectorSpace(lhs_h_->NRows());
  SmartPtr<DenseVector> retval_x_nonconst = retval_nc_space->MakeNewDenseVector();
  lhs_h_->MultVector(1.0,*retval_x,0.0,*retval_x_nonconst);
  retval_nc_space = new DenseVectorSpace(lhs_c_->NRows());
  SmartPtr<DenseVector> retval_c_nonconst = retval_nc_space->MakeNewDenseVector();
  lhs_c_->MultVector(1.0,*retval_x,0.0,*retval_c_nonconst);
  retval_nc_space = new DenseVectorSpace(lhs_d_->NRows());
  SmartPtr<DenseVector> retval_d_nonconst = retval_nc_space->MakeNewDenseVector();
  lhs_d_->MultVector(1.0,*retval_x,0.0,*retval_d_nonconst);

  // reduce to the desired quantities
  retval_x = shrink(retval_x_nonconst,shift_x_indices_);
  retval_c = shrink(retval_c_nonconst,shift_c_indices_);
  retval_d = shrink(retval_d_nonconst,shift_d_indices_);

  //setup return value
  VectorSet retval;
  retval.top = target.top->MakeNewIteratesVectorCopy();
  retval.x = retval_x;
  retval.y_c = retval_c;
  retval.y_d = retval_d;

  return retval;
}

VectorSet LinKKTFaster::undoConditionning(SmartPtr<IpoptApplication> app, const VectorSet& target) const
{
  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(target.top->Dim()+target.x->Dim()+target.y_c->Dim()+target.y_d->Dim());

  VectorSet retval;
  retval.top = target.top->MakeNewIteratesVectorCopy();
  retval.x = dynamic_cast<DenseVector*>(target.x->MakeNewCopy());
  retval.y_c = dynamic_cast<DenseVector*>(target.y_c->MakeNewCopy());
  retval.y_d = dynamic_cast<DenseVector*>(target.y_d->MakeNewCopy());

  SmartPtr<IteratesVector> retval_top = retval.top->MakeNewIteratesVectorCopy();

  // get KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

  // do actual backsolve
  pd_solver->Solve(1.0, 0.0, *target.top, *retval_top);

  SmartPtr<const Vector> retval_x = retval_top->x();
  SmartPtr<const Vector> retval_c = retval_top->y_c();
  SmartPtr<const Vector> retval_d = retval_top->y_d();

  retval.top->Set(0.0);
  retval.top->Set_x(*retval_x);
  retval.top->Set_y_c(*retval_c);
  retval.top->Set_y_d(*retval_d);

 return retval;
}

SmartPtr<const DenseVector> LinKKTFaster::computeDeltaP(const IntervalInfoSet& intervals,const Index& interval) const
{
  const Index n_parameters = intervals.getParameterCount();
  Index counter = 0;
  std::vector<Index> parameterflags = intervals.getParameterIDVec();
  std::vector<Index> intervalflags = intervals.getIntervalIDVec();
  std::vector<Index> parameterIDs = intervals.getParameterIDs();
  IntWidthScaling *int_widthscaler = new IntWidthScaling();
  std::vector<Number> int_widths = int_widthscaler->scaleIntervalWidths(intervals,interval);
  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(2*n_parameters);
  SmartPtr<DenseVector> retval = retval_space->MakeNewDenseVector();
  Number* ret_values = new Number[2*n_parameters];

  for (int k=0; k<intervalflags.size(); k++) {
    for (int j=0; j<parameterIDs.size(); j++) {
      if (intervalflags[k] == interval) {
	if (parameterflags[k] == parameterIDs[j]) {
	  if (intervals.isUpper(k))
	    ret_values[counter] = -0.5*int_widths[j];
	   else
	    ret_values[counter] = 0.5*int_widths[j];
	  counter++;
	}
      }
    }
  }

  retval->SetValues(ret_values);
  return dynamic_cast<const DenseVector*>(GetRawPtr(retval));
}

SplitDecision LinKKTFaster::decideInterval(const std::vector<SplitApproximation>& approximates) const
{
  /* waiting to be implemeted*/
}

SplitApproximation LinKKTFaster::applyAlgorithmOnInterval(SmartPtr<IpoptApplication> app, const Index& interval)
{
  // initialize neccessary quantities
  const Index n_parameters = intervals_.getParameterCount();
  std::vector<Index> intervalIDs = intervals_.getIntervalIDVec();
  SmartPtr<DenseVectorSpace> z_space = new DenseVectorSpace(rhs_dim_);
  SmartPtr<DenseVector> z_postsplit = z_space->MakeNewDenseVector();
  SmartPtr<const DenseVector> deltap;

  SmartPtr<MultiVectorMatrixSpace> sense_space = new MultiVectorMatrixSpace(n_parameters*2,*z_space);
  SmartPtr<MultiVectorMatrix> sense = sense_space->MakeNewMultiVectorMatrix();
  SmartPtr<DenseVectorSpace> res_space = new DenseVectorSpace(sense->NRows());
  SmartPtr<DenseVector> res_sense = res_space->MakeNewDenseVector();
  VectorSet rhs;
  VectorSet rhs_scaled;
  VectorSet lhs;

  // start algorithm on each existing intervalID
  Index s_count = 0;
  for (Index j=0;j<intervalIDs.size();j++) {
    if (intervalIDs[j]==interval) {

  // assign private attributes depending on chosen interval
      if (assignIntAndParaDepAttr(interval,j)) {
	// create VectorSet for computation calls

	SmartPtr<IteratesVector> rhs_top = app->IpoptDataObject()->curr()->MakeNewIteratesVector();
	rhs_top->Set(0.0);
	rhs_top->Set_x(*x_i_);
	rhs_top->Set_y_c(*y_c_i_);
	rhs_top->Set_y_d(*y_d_i_);
	rhs.top = rhs_top;
	rhs.x = dynamic_cast<DenseVector*>(rhs_shift_h_->MakeNewCopy());
	rhs.y_c = dynamic_cast<DenseVector*>(rhs_shift_c_->MakeNewCopy());
	rhs.y_d = dynamic_cast<DenseVector*>(rhs_shift_d_->MakeNewCopy());
	// printVectorSet(app,rhs,"rhs");
	z_postsplit = transformVectorSet(computeS(app,rhs));
	z_postsplit->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"S_fast");
	// compute solution of the GMRES-Step
	const Number steplength = computeRMultS(app,rhs)*computeSMultS(app,rhs);
	/*      printf("\n\n");
		char buffer[16];
		sprintf(buffer,"rhs_i_[%d]",j);
		rhs_i_->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, buffer);
	*/

	//  printf("\n n_u_= %d, n_x_ = %d, n_c_ = %d, n_d_ = %d", n_u_, n_x_, n_c_, n_d_);
	//	printVectorSet(app,rhs,"rhs");
	if (steplength) {
	  //	printf("\nLinKKTFaster::applyAlgorithmOnInterval(): MINRES-Step is %e.",steplength);
	  rhs_scaled = scaleVectorSet(rhs,steplength);
	  // z_postsplit = transformVectorSet(rhs_scaled);
	  // z_postsplit->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"z_sc_fast");
	  // z_postsplit = transformVectorSet(computeS(app,rhs));
	  // z_postsplit->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"S_fast");
	//printVectorSet(app,rhs,"rhs");
//	printVectorSet(app,rhs_scaled,"rhs_scaled");
	} else
	  printf("\nLinKKTFaster::applyAlgorithmOnInterval(): ERROR: MINRES-Step is 0.");

	// undo conditionning step to get the approximation of the postsplit sensitivity column
	lhs = undoConditionning(app,rhs_scaled);
//	printVectorSet(app,rhs_scaled,"rhs_scaled");
//	printVectorSet(app,lhs,"lhs");
	// compare sense vectors for test purposes
	//      SmartPtr<Vector> compy = z_postsplit->MakeNew();
	/*      printf("\n");
		for (int i=0;i<n_x_+n_c_+n_d_+n_u_;i++)
		printf("\n different");
		//      compy->AddTwoVectors(1.0,*z_copy,-1.0,*z_postsplit_cond,0.0);
		printf("\n\n");
		z_postsplit->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"z_post");
		//      compy->AddTwoVectors(1.0,*z_copy,-1.0,*z_postsplit,0.0);
		printf("\n\n");
		//      compy->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"compycomp");

		*/

	// transform the VectorSet into a DenseVector
//	printVectorSet(app,lhs,"lhs");
	z_postsplit = transformVectorSet(lhs);
//	printVectorSet(app,lhs,"lhs");
//	z_postsplit->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"z_post_f");

	// add the resulting sense vector to matrix ....
	sense->SetVector(s_count,*z_postsplit);
	s_count++;

      } else
	printf("\nLinKKTFaster::applyAlgorithmOnInterval(): ERROR: Unable to assign interval dependent attributes!");
    }
  }
  // get the parameter perturbation for this split
  deltap = computeDeltaP(intervals_,interval);

  // calculate absolute sensitivity values
  sense->MultVector(1.0,*deltap,0.0,*res_sense);

  printf("\n\n");
  res_sense->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"result_f");
  printf("\n\n");
  //  std::vector<Index> x_intervalIDs = (dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_->OwnerSpace())))->GetIntegerMetaData("intervalID");
  //  printf("\n");
  // for (int i=0;i<x_intervalIDs.size();i++)
  //   printf("\n interval---%d----",x_intervalIDs[i]);


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

void LinKKTFaster::printVectorSet(SmartPtr<IpoptApplication>app,const VectorSet& target, const std::string& label) const
{
  char buffer[63];
  sprintf(buffer,"%s.top.x",label.c_str());
  printf("\n");
  target.top->x()->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, buffer);
  printf("\n");
  sprintf(buffer,"%s.top.y_c",label.c_str());
  target.top->y_c()->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, buffer);
  printf("\n");
  sprintf(buffer,"%s.top.y_d",label.c_str());
  target.top->y_d()->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, buffer);
  printf("\n");
  sprintf(buffer,"%s.x",label.c_str());
  target.x->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, buffer);
  sprintf(buffer,"%s.y_c",label.c_str());
  target.y_c->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, buffer);
  sprintf(buffer,"%s.y_d",label.c_str());
  target.y_d->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, buffer);
}

VectorSet LinKKTFaster::scaleVectorSet(const VectorSet& target,const Number& factor) const
{
  VectorSet retval;
  retval.top = target.top->MakeNewIteratesVectorCopy();
  retval.x = dynamic_cast<DenseVector*>(target.x->MakeNewCopy());
  retval.y_c = dynamic_cast<DenseVector*>(target.y_c->MakeNewCopy());
  retval.y_d = dynamic_cast<DenseVector*>(target.y_d->MakeNewCopy());

  retval.top->Scal(factor);
  retval.x->Scal(factor);
  retval.y_c->Scal(factor);
  retval.y_d->Scal(factor);

  return retval;
}

SmartPtr<DenseVector> LinKKTFaster::transformVectorSet(const VectorSet& target) const
{
  // initialize
  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(rhs_dim_);
  SmartPtr<DenseVector> retval = retval_space->MakeNewDenseVector();
  Index abs_pos = 0;
  SmartPtr<DenseVector> top_part = retval->MakeNewDenseVector();
  SmartPtr<DenseVector> x_part = retval->MakeNewDenseVector();
  SmartPtr<DenseVector> y_c_part = retval->MakeNewDenseVector();
  SmartPtr<DenseVector> y_d_part = retval->MakeNewDenseVector();

  // expand parts
  top_part = expand(target.top,retval_space->Dim(),abs_pos);
  abs_pos+= n_x_+n_c_+n_u_+n_d_;
  x_part = expand(target.x,retval->Dim(),abs_pos);
  abs_pos+= target.x->Dim();
  y_c_part = expand(target.y_c,retval->Dim(),abs_pos);
  abs_pos+= target.y_c->Dim();
  y_d_part = expand(target.y_d,retval->Dim(),abs_pos);

  // sum parts
  retval->AddTwoVectors(1.0,*top_part,1.0,*x_part,0.0);
  retval->AddTwoVectors(1.0,*y_d_part,1.0,*y_c_part,1.0);

  return retval;
}

/** given an Index-type vector containing a list of intervalIDs, loads indices matching the given interval into the shift_indices part, the rest into static_indices **/
bool LinKKTFaster::splitIntervalIndices(const std::vector<Index>& intervalIDs,std::vector<Index>& static_indices,std::vector<Index>& shift_indices,const Index& interval)
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

  //  SplitAlgorithm* splitter = assignSplitAlgorithm(app);
  SplitAlgorithm* splitter = new LinKKTFaster(app);
  SplitDecision resulting_split = splitter->applySplitAlgorithm(app);
  // printf("\nCHANGING SPLITMODES.");
  // splitter = new LinearizeKKTWithMINRES(app);
  // resulting_split = splitter->applySplitAlgorithm(app);

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

