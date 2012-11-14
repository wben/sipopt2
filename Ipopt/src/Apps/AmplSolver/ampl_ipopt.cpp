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

/* forward declaration*/
class ShiftVector;
class LinearizeKKT : public SplitAlgorithm
{
public:
  LinearizeKKT(SmartPtr<IpoptApplication> app);
  /* structural parts*/
  SplitDecision applySplitAlgorithm(SmartPtr<IpoptApplication> app);
  SplitApproximation applyAlgorithmOnInterval(SmartPtr<IpoptApplication> app, const Index& interval);
  SplitDecision decideInterval(const std::vector<SplitApproximation>& approximates) const;

  /*handling and manipulating data specifically*/
  bool splitIntervalIndices(const std::vector<Index>& intervalIDs,std::vector<Index>& shift_indices,const Index& interval);
  SmartPtr<DenseVector> extractColumn(SmartPtr<const Matrix> original,const Index& column) const;
  SmartPtr<DenseVector> expand(SmartPtr<DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const;
  SmartPtr<DenseVector> expand(SmartPtr<DenseVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<DenseVector> expand(SmartPtr<IteratesVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<DenseVector> shrink(SmartPtr<DenseVector> original,const std::vector<Index>& indices) const;
  /*algorithm specific components*/
  Number computeRMultS(SmartPtr<IpoptApplication> app,SmartPtr<ShiftVector> target) const;
  Number computeSMultS(SmartPtr<IpoptApplication> app,SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeS(SmartPtr<IpoptApplication> app,SmartPtr<ShiftVector> target) const;
  /*  SmartPtr<const DenseVector> computeR() const;*/ // not neccessary for x_0=0 and 0 steps.
  SmartPtr<ShiftVector> computeDMultYrhs(SmartPtr<IpoptApplication> app,SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> undoConditionning(SmartPtr<IpoptApplication> app, SmartPtr<ShiftVector> target) const;
  SmartPtr<const DenseVector> computeDeltaP(const IntervalInfoSet& intervals,const Index& interval) const;

private:
  //original Ipopt data
  SmartPtr<const Matrix> Wp_;
  SmartPtr<const Matrix> Ap_;
  SmartPtr<const Matrix> Bp_;
  SmartPtr<const SymMatrix> W_;
  SmartPtr<const Matrix> A_;
  SmartPtr<const Matrix> B_;

  // extracted data - not interval specific
  IntervalInfoSet intervals_;
  std::vector<Index> u_indices_;
  std::vector<Index> x_intervalIDs_;
  std::vector<Index> c_intervalIDs_;
  std::vector<Index> d_intervalIDs_;
  Index n_i_;
  Index n_u_;
  Index rhs_dim_;
  Index top_dim_;

  // extracted data - interval specific
  SmartPtr<DenseVector> u_i_;
  std::vector<Index> shift_x_indices_;
  std::vector<Index> shift_c_indices_;
  std::vector<Index> shift_d_indices_;

  // manipulated or constructed data - interval specific
  SmartPtr<DenseVector> x_i_;
  SmartPtr<DenseVector> y_c_i_;
  SmartPtr<DenseVector> y_d_i_;
};

/*forward declaration*/
//class ShiftVectorSpace;

class ShiftVector : public ReferencedObject
{    /*NOTICE: at the moment, all parts of the implementation that resemble Vector base class
implementation will take action on the local DenseVector representation of the shift Vector - hence
work on the top->x(), top->y_c(), top->y_d(), x(),y_c() and y_d() part of the ShiftVector. No other
parts of the local IteratesVector are implied. */
public:
  ShiftVector(SmartPtr<IteratesVector> top,
	      SmartPtr<DenseVector> x,
	      SmartPtr<DenseVector>y_c,
	      SmartPtr<DenseVector> y_d);
  ShiftVector(const ShiftVector& rhs);
  ShiftVector& operator=(const ShiftVector& rhs);
  SmartPtr<DenseVector> getDVector() const;
  void Scal(const Number& factor);
  //  SmartPtr<ShiftVectorSpace> OwnerSpace() const;
  void Print(SmartPtr<const Journalist> jnlst,
               EJournalLevel level,
               EJournalCategory category,
               const std::string& name,
               Index indent=0,
               const std::string& prefix="") const;
  void Print(const Journalist& jnlst,
               EJournalLevel level,
               EJournalCategory category,
               const std::string& name,
               Index indent=0,
               const std::string& prefix="") const;
  Number Dot(const Vector &x) const;
  Number Dot(const ShiftVector &x) const;
  Index Dim() const;
  void Set(Number alpha);
  void AddOneVector(const Number& a,const ShiftVector& v1,const Number& c);
  void AddTwoVectors(const Number& a,const ShiftVector& v1,const Number& b,const ShiftVector& v2,const Number& c);
  SmartPtr<ShiftVector> MakeNewShiftVector() const;
  //  SmartPtr<Vector> MakeNew() const; // disabled: no Vector child yet
  SmartPtr<IteratesVector>top() const;
  SmartPtr<DenseVector>x() const;
  SmartPtr<DenseVector>y_c() const;
  SmartPtr<DenseVector>y_d() const;

private:
  ShiftVector();
  void top(SmartPtr<IteratesVector>top); //remember to set privates
  void x(SmartPtr<DenseVector>x);  //remember to set privates
  void y_c(SmartPtr<DenseVector>y_c);  //remember to set privates
  void y_d(SmartPtr<DenseVector>y_d);  //remember to set privates
  //  void ownerspace(SmartPtr<ShiftVectorSpace>ospace);  //remember to set privates
  void exp(SmartPtr<DenseVector> exp);    //remember to set privates + space

  SmartPtr<IteratesVector> top_;
  SmartPtr<DenseVector> x_;
  SmartPtr<DenseVector> y_c_;
  SmartPtr<DenseVector> y_d_;
  SmartPtr<DenseVectorSpace> exp_space_;
  SmartPtr<DenseVector> exp_;
  //  SmartPtr<ShiftVectorSpace> owner_space_;
  Index dim_;
  Index top_dim_;
  Index x_dim_;
  Index y_c_dim_;
  Index y_d_dim_;

};
/*
class ShiftVectorSpace : public ReferencedObject
{
public:
  ShiftVectorSpace(SmartPtr<IteratesVectorSpace> top,
		   SmartPtr<DenseVectorSpace> x,
		   SmartPtr<DenseVectorSpace>y_c,
		   SmartPtr<DenseVectorSpace> y_d);
  //  SmartPtr<Vector> MakeNew() const; // disabled: no Vector child yet
  SmartPtr<ShiftVector> MakeNewShiftVector() const;
  Index Dim() const;
  Index top_Dim() const;
  Index x_Dim() const;
  Index y_c_Dim() const;
  Index y_d_Dim() const;
private:
  SmartPtr<IteratesVectorSpace> top_space_;
  SmartPtr<DenseVectorSpace> x_space_;
  SmartPtr<DenseVectorSpace> y_c_space_;
  SmartPtr<DenseVectorSpace> y_d_space_;
  Index dim_;
  Index top_dim_;
  Index x_dim_;
  Index y_c_dim_;
  Index y_d_dim_;

};
*/
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

ShiftVector::ShiftVector(SmartPtr<IteratesVector> top, SmartPtr<DenseVector> x, SmartPtr<DenseVector>y_c, SmartPtr<DenseVector> y_d)
{
  assert((GetRawPtr(top)) && (GetRawPtr(x)) && (GetRawPtr(y_c)) && (GetRawPtr(y_d)));

  top_ = top;
  x_ = x;
  y_c_ = y_c;
  y_d_ = y_d;


  top_dim_ = top_->x()->Dim()+top_->y_c()->Dim()+top_->y_d()->Dim();
  x_dim_ = x_->Dim();
  y_c_dim_ = y_c_->Dim();
  y_d_dim_ = y_d_->Dim();
  dim_ = top_dim_+x_dim_+y_c_dim_+y_d_dim_;

  //  owner_space_ = new ShiftVectorSpace(tmp_space,x_->OwnerSpace(),y_c_->OwnerSpace(),y_d_->OwnerSpace());

  // construct DenseVector expansion of the stored values
  exp_space_ = new DenseVectorSpace(dim_);
  exp_ = exp_space_->MakeNewDenseVector();

  Number* exp_values = new Number[dim_];
  Index abs_pos = 0;

  const Number* top_x_val = dynamic_cast<const DenseVector*>(GetRawPtr(top->x()))->Values();
  const Number* top_y_c_val = dynamic_cast<const DenseVector*>(GetRawPtr(top->y_c()))->Values();
  const Number* top_y_d_val = dynamic_cast<const DenseVector*>(GetRawPtr(top->y_d()))->Values();

  const Number* x_val = x_->Values();
  const Number* y_c_val = y_c_->Values();
  const Number* y_d_val = y_d_->Values();

  const Index top_x_dim = top_->x()->Dim();
  const Index top_y_c_dim = top_->y_c()->Dim();
  const Index top_y_d_dim = top_->y_d()->Dim();

  for (int i=0;i<top_x_dim;i++){
    exp_values[abs_pos] = 0;
    if (top_x_val)
      exp_values[abs_pos] = top_x_val[i];
    abs_pos++;
  }
  for (int i=0;i<top_y_c_dim;i++){
    exp_values[abs_pos] = 0;
    if (top_y_c_val)
      exp_values[abs_pos] = top_y_c_val[i];
    abs_pos++;
  }
  for (int i=0;i<top_y_d_dim;i++){
    exp_values[abs_pos] = 0;
    if (top_y_d_val)
      exp_values[abs_pos] = top_y_d_val[i];
    abs_pos++;
  }
  for (int i=0;i<x_dim_;i++){
    exp_values[abs_pos] = 0;
    if (x_val)
      exp_values[abs_pos] = x_val[i];
    abs_pos++;
  }
  for (int i=0;i<y_c_dim_;i++){
    exp_values[abs_pos] = 0;
    if (y_c_val)
      exp_values[abs_pos] = y_c_val[i];
    abs_pos++;
  }
  for (int i=0;i<y_d_dim_;i++){
    exp_values[abs_pos] = 0;
    if (y_d_val)
      exp_values[abs_pos] = y_d_val[i];
    abs_pos++;
  }
  exp_->SetValues(exp_values);
}

ShiftVector::ShiftVector(const ShiftVector& rhs)
{
  top_ = rhs.top();
  x_ = rhs.x();
  y_c_ = rhs.y_c();
  y_d_ = rhs.y_d();

  exp_ = rhs.getDVector();
  exp_space_ = new DenseVectorSpace(exp_->Dim());

  //  owner_space_ = new ShiftVectorSpace(top_->OwnerSpace,x_->OwnerSpace(),y_c_->OwnerSpace(),y_d_->OwnerSpace());

  dim_ = exp_->Dim();
  top_dim_ = top_->x()->Dim()+top_->y_c()->Dim()+top_->y_d()->Dim();
  x_dim_ = x_->Dim();
  y_c_dim_ = y_c_->Dim();
  y_d_dim_ = y_d_->Dim();
}

ShiftVector& ShiftVector::operator=(const ShiftVector &rhs)
{
  if (this!= &rhs)
    {
      this->top(rhs.top());
      this->x(rhs.x());
      this->y_c(rhs.y_c());
      this->y_d(rhs.y_d());
      //  this->ownerspace(rhs.OwnerSpace());
      this->exp(rhs.getDVector());
    }
  return *this;
}

SmartPtr<DenseVector> ShiftVector::getDVector() const
{
  return exp_;
}

void ShiftVector::Scal(const Number& factor)
{
  SmartPtr<DenseVector> top_x = dynamic_cast<DenseVector*>(top_->x()->MakeNewCopy());
  SmartPtr<DenseVector> top_y_c = dynamic_cast<DenseVector*>(top_->y_c()->MakeNewCopy());
  SmartPtr<DenseVector> top_y_d = dynamic_cast<DenseVector*>(top_->y_d()->MakeNewCopy());

  top_x->Scal(factor);
  top_y_c->Scal(factor);
  top_y_d->Scal(factor);

  top_->Set_x_NonConst(*top_x);
  top_->Set_y_c_NonConst(*top_y_c);
  top_->Set_y_d_NonConst(*top_y_d);

  x_->Scal(factor);
  y_c_->Scal(factor);
  y_d_->Scal(factor);
  exp_->Scal(factor);
}
/*
SmartPtr<ShiftVectorSpace> ShiftVector::OwnerSpace() const
{
  return owner_space_;
  } */

void ShiftVector::Print(SmartPtr<const Journalist> jnlst, EJournalLevel level, EJournalCategory category, const std::string& name, Index indent, const std::string& prefix) const
{
  exp_->Print(jnlst,level,category,name,indent,prefix);
}

void ShiftVector::Print(const Journalist& jnlst, EJournalLevel level, EJournalCategory category, const std::string& name, Index indent, const std::string& prefix) const
{
  exp_->Print(jnlst,level,category,name,indent,prefix);
}

Number ShiftVector::Dot(const Vector &x) const
{
  assert(x.Dim()==dim_);
  return exp_->Dot(x);
}

Number ShiftVector::Dot(const ShiftVector &x) const
{
  assert(x.Dim()==dim_);
  SmartPtr<DenseVector> exp = x.getDVector();
  return exp_->Dot(*exp);
}
void ShiftVector::AddOneVector(const Number& a,const ShiftVector& v1,const Number& c)
{
  AddTwoVectors(a,v1,0.0,*this,c);
}

void ShiftVector::AddTwoVectors(const Number& a,const ShiftVector& v1,const Number& b,const ShiftVector& v2,const Number& c)
{
  if (!c)
    Set(0.0);

  if (a && b) {
    SmartPtr<const Vector> top1_x = v1.top()->x();
    SmartPtr<const Vector> top1_y_c = v1.top()->y_c();
    SmartPtr<const Vector> top1_y_d = v1.top()->y_d();

    SmartPtr<DenseVector> x1 = v1.x();
    SmartPtr<DenseVector> y_c1 = v1.y_c();
    SmartPtr<DenseVector> y_d1 = v1.y_d();

    SmartPtr<DenseVector> x2 = v2.x();
    SmartPtr<DenseVector> y_c2 = v2.y_c();
    SmartPtr<DenseVector> y_d2 = v2.y_d();

    SmartPtr<const Vector> top2_x = v2.top()->x();
    SmartPtr<const Vector> top2_y_c = v2.top()->y_c();
    SmartPtr<const Vector> top2_y_d = v2.top()->y_d();

    SmartPtr<Vector> x = top()->x()->MakeNewCopy();
    SmartPtr<Vector> y_c = top()->y_c()->MakeNewCopy();
    SmartPtr<Vector> y_d = top()->y_d()->MakeNewCopy();

    x->AddTwoVectors(a,*top1_x,b,*top2_x,c);
    y_c->AddTwoVectors(a,*top1_y_c,b,*top2_y_c,c);
    y_d->AddTwoVectors(a,*top1_y_d,b,*top2_y_d,c);

    top_->Set_x_NonConst(*x);
    top_->Set_y_c_NonConst(*y_c);
    top_->Set_y_d_NonConst(*y_d);

    x_->AddTwoVectors(a,*x1,b,*x2,c);
    y_c_->AddTwoVectors(a,*y_c1,b,*y_c2,c);
    y_d_->AddTwoVectors(a,*y_d1,b,*y_d2,c);
  }

}

Index ShiftVector::Dim() const
{
  return dim_;
}

void ShiftVector::Set(Number alpha)
{
  top_->Set(alpha);
  x_->Set(alpha);
  y_c_->Set(alpha);
  y_d_->Set(alpha);
  exp_->Set(alpha);
}

SmartPtr<ShiftVector> ShiftVector::MakeNewShiftVector() const
{
  SmartPtr<IteratesVector> top = dynamic_cast<IteratesVector*>(top_->MakeNew());;
  SmartPtr<DenseVector> x = x_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_c = y_c->MakeNewDenseVector();
  SmartPtr<DenseVector> y_d = y_d->MakeNewDenseVector();
  return  new ShiftVector(top,x,y_c,y_d);
}
/*  // disabled: no Vector child yet
SmartPtr<Vector> ShiftVector::MakeNew() const
{

  return exp_->MakeNew();
  }*/

SmartPtr<IteratesVector> ShiftVector::top() const
{
  return top_;
}

SmartPtr<DenseVector> ShiftVector::x() const
{
  return x_;
}

SmartPtr<DenseVector> ShiftVector::y_c() const
{
  return y_c_;
}

SmartPtr<DenseVector> ShiftVector::y_d() const
{
  return y_d_;
}

void ShiftVector::top(SmartPtr<IteratesVector>top)
{
  top_dim_ = top->Dim();
  top_ = top;
}

void ShiftVector::x(SmartPtr<DenseVector>x)
{
  x_dim_ = x->Dim();
  x_ = x;
}

void ShiftVector::y_c(SmartPtr<DenseVector>y_c)
{
  y_c_dim_ = y_c->Dim();
  y_c_ = y_c;
}

void ShiftVector::y_d(SmartPtr<DenseVector>y_d)
{
  y_d_dim_ = y_d->Dim();
  y_d_ = y_d;
}
/*
void ShiftVector::ownerspace(SmartPtr<ShiftVectorSpace>ospace)
{
  owner_space_ = ospace;
  dim_ = ospace->Dim();
  }*/

void ShiftVector::exp(SmartPtr<DenseVector> exp)
{
  exp_ = exp;
  exp_space_ = new DenseVectorSpace(exp_->Dim());
}
/*
ShiftVectorSpace::ShiftVectorSpace(SmartPtr<IteratesVectorSpace> top,SmartPtr<DenseVectorSpace> x,SmartPtr<DenseVectorSpace>y_c, SmartPtr<DenseVectorSpace> y_d)
{
  assert(GetRawPtr(top) && GetRawPtr(x) && GetRawPtr(y_c) && GetRawPtr(y_d));

  //  SmartPtr<IteratesVectorSpace>  tmp_space = dynamic_cast<IteratesVectorSpace*>(GetRawPtr(top));

  top_space_ = top;
  x_space_ = x;
  y_c_space_ = y_c;
  y_d_space_= y_d;

  top_dim_ = x_space_->Dim()+y_c_space_->Dim()+y_d_space_->Dim();
  x_dim_ = x->Dim();
  y_c_dim_ = y_c->Dim();
  y_d_dim_ = y_d->Dim();
  dim_ = top_dim_+x_dim_+y_c_dim_+y_d_dim_;
}

/*  // disabled: no Vector child yet
SmartPtr<Vector>ShiftVectorSpace::MakeNew() const
{
  SmartPtr<IteratesVector> top = top_space_->MakeNewIteratesVector();
  SmartPtr<DenseVector> x = x_space_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_c = y_c_space_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_d = y_d_space_->MakeNewDenseVector();

  return dynamic_cast<Vector*>(new ShiftVector(top,x,y_c,y_d));
}   */
 /*
SmartPtr<ShiftVector> ShiftVectorSpace::MakeNewShiftVector() const
{
  SmartPtr<IteratesVector> top = top_space_->MakeNewIteratesVector();
  SmartPtr<DenseVector> x = x_space_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_c = y_c_space_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_d = y_d_space_->MakeNewDenseVector();

  return new ShiftVector(top,x,y_c,y_d);
}

Index ShiftVectorSpace::Dim() const
{
  return dim_;
}

Index ShiftVectorSpace::top_Dim() const
{
  return top_dim_;
}

Index ShiftVectorSpace::x_Dim() const
{
  return x_dim_;
}

Index ShiftVectorSpace::y_c_Dim() const
{
  return y_c_dim_;
}

Index ShiftVectorSpace::y_d_Dim() const
{
  return y_d_dim_;
}
 */

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

/*SmartPtr<MultiVectorMatrix> getMVMFromMatrix(SmartPtr<const Matrix> matrix)
{
  // set up neccessary Spaces and return value candidates
  const Index n_cols = (matrix->OwnerSpace())->NCols();
  const Index n_rows = (matrix->OwnerSpace())->NRows();

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
  } */


SplitDecision SplitWRTControlSensitivities::applySplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  SmartPtr<OptionsList> options = app->Options();
  BranchingCriterion* branchmode = assignBranchingMethod(options);
  std::vector<SplitChoice> splitchoices = branchmode->branchSensitivityMatrix(app);

  ControlSelector* pickfirst = assignControlMethod(options);
  return pickfirst->decideSplitControl(splitchoices);
}

SplitAlgorithm* assignSplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  SmartPtr<OptionsList> options = app->Options();
  // read from options what kind of controlselector to return
  std::string sensemode;
  SplitAlgorithm* retval = new LinearizeKKT(app);
  if (options->GetStringValue("sensemode",sensemode,"")) {
    if (sensemode == "control") {
      retval = new SplitWRTControlSensitivities();
    } if (sensemode == "GMRES") {
      retval = new LinearizeKKT(app);
    }
  } else
    printf("\nassignSplitAlgorithm(): ERROR: No splitAlgorithm chosen!");
  return retval;
}

LinearizeKKT::LinearizeKKT(SmartPtr<IpoptApplication> app)
{
  x_intervalIDs_.clear();
  c_intervalIDs_.clear();
  d_intervalIDs_.clear();

  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  //assign local original Ipopt data
  SmartPtr<const DenseVector> x_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->x()));
SmartPtr<const DenseVector>  p_ = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
SmartPtr<const DenseVector>  y_c_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->y_c()));
SmartPtr<const DenseVector>  y_d_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->y_d()));
  Wp_ = orig_nlp->h_p(*x_, 1.0, *y_c_, *y_d_);
  Ap_  = orig_nlp->jac_c_p(*x_);
  Bp_  = orig_nlp->jac_d_p(*x_);
  W_ = orig_nlp->h(*x_, 1.0, *y_c_, *y_d_);
  A_ = orig_nlp->jac_c(*x_);
  B_ = orig_nlp->jac_d(*x_);

  SmartPtr<const DenseVectorSpace> x_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x_->OwnerSpace()));
  SmartPtr<const DenseVectorSpace> c_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(y_c_->OwnerSpace()));
  SmartPtr<const DenseVectorSpace> d_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(y_d_->OwnerSpace()));

  //assign local extracted data - non interval specific
  intervals_ = IntervalInfoSet(p_);
  std::vector<Index> x_intervalIDs = x_space->GetIntegerMetaData("intervalID");
  std::vector<Index> u_indices;
  std::vector<Index> x_indices;
  if (splitIntervalIndices(x_intervalIDs,u_indices,0)) {
    assert(u_indices.size());
    u_indices_ = u_indices;
    x_intervalIDs_=x_intervalIDs;
    }

  if (c_space->HasIntegerMetaData("intervalID")){
    c_intervalIDs_ = c_space->GetIntegerMetaData("intervalID");
  } else
    printf("LinearizeKKT::LinearizeKKT(): Error - no intervalIDs for c_space!\n");
  if (d_space->HasIntegerMetaData("intervalID")) {
    d_intervalIDs_ = d_space->GetIntegerMetaData("intervalID");
  } else
    printf("LinearizeKKT::LinearizeKKT(): Error - no intervalIDs for d_space!\n");

  n_i_ = intervals_.getIntervalCount();
  n_u_ = u_indices_.size();
  Index  n_x_ = x_intervalIDs_.size()-n_u_;
  Index  n_c_ = c_intervalIDs_.size();
  Index  n_d_ = d_intervalIDs_.size();

  top_dim_ = n_x_+n_u_+n_c_+n_d_;
  rhs_dim_ = top_dim_+int((n_x_+n_c_+n_d_)/n_i_);
}

/* apply GMRES-approximation of a split to all intervals and decide for the smartest split*/
SplitDecision LinearizeKKT::applySplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  std::vector<SplitApproximation> approximates(n_i_);
  for (int i=0;i<n_i_;i++) {

      // get GMRES approximated split results wrt each single interval
      // i+1: in case intervalIDs start with 1 (which they do)
      approximates[i] = applyAlgorithmOnInterval(app,i+1);
  }

  SplitDecision retval;

  // chose the split with best results
  retval = this->decideInterval(approximates);

  ///////////////////////only for the sake of a working python/ampl interface////
  ////////manually performing a ctrlwise branch after LinKKT/////////////////////
  SmartPtr<OptionsList> options = app->Options();
  BranchingCriterion* branchmode = assignBranchingMethod(options);
  std::vector<SplitChoice> splitchoices = branchmode->branchSensitivityMatrix(app);
  ControlSelector* pickfirst = assignControlMethod(options);
  retval = pickfirst->decideSplitControl(splitchoices);

  return retval;
}

// get a specific column of a given const Matrix* as a const Densevector*
SmartPtr<DenseVector> LinearizeKKT::extractColumn(SmartPtr<const Matrix> original,const Index& column) const
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
  SmartPtr<DenseVector> retval = retval_space->MakeNewDenseVector();
  original->MultVector(1.0,*unit,0.0,*retval);

  return retval;
}

SmartPtr<DenseVector> LinearizeKKT::expand(SmartPtr<DenseVector> original,const Index& large_dim, const Index& start_idx) const
{
  const Index small_dim = original->Dim();
  assert(large_dim>=small_dim);

  Index* exppos = new Index[small_dim];
  for (int i=0;i<small_dim;i++)
    exppos[i]=i+start_idx;

  SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
  SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
  SmartPtr<DenseVectorSpace> expanded_space = new DenseVectorSpace(large_dim);
  SmartPtr<DenseVector> new_vec = expanded_space->MakeNewDenseVector();

  em->MultVector(1.0,*original,0.0,*new_vec);

  return new_vec;
}

SmartPtr<DenseVector> LinearizeKKT::expand(SmartPtr<IteratesVector> original,const Index& large_dim, const Index& start_idx) const
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
SmartPtr<DenseVector> LinearizeKKT::expand(SmartPtr<DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const
{
  const Index small_dim = indices.size();

  assert(large_dim>=small_dim);

  if (large_dim==small_dim) {
//    printf("\nLinearizeKKT::expandVector(): WARNING: called to expand a vector to original size.");
    return original;
  } else {
    Index* exppos = new Index[small_dim];
    for (int i=0;i<small_dim;i++)
      exppos[i]=indices[i];

    SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
    SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
    SmartPtr<DenseVectorSpace> expanded_space = new DenseVectorSpace(large_dim);
    SmartPtr<DenseVector> new_vec = expanded_space->MakeNewDenseVector();

    em->MultVector(1.0,*original,0.0,*new_vec);

    return new_vec;
  }
  printf("\nLinearizeKKT::expand(): Unknown ERROR.");
  return NULL;
}

/* shrink Vector original from larger (old) dim to smaller dim, with the vector indices listing which indices from the old vector to keep*/
SmartPtr<DenseVector> LinearizeKKT::shrink(SmartPtr<DenseVector> original,const std::vector<Index>& indices) const
{
  const Index large_dim = original->Dim();
  const Index small_dim = indices.size();

  assert(large_dim>=small_dim);
  if (large_dim==small_dim) {
    //    printf("\nLinearizeKKT::shrink(): WARNING: called to shrink a vector to original size.");
    return original;
  } else {

    Index* exppos = new Index[small_dim];
    for (int i=0;i<small_dim;i++)
      exppos[i]=indices[i];


    SmartPtr<ExpansionMatrixSpace> em_space = new ExpansionMatrixSpace(large_dim,small_dim,exppos);
    SmartPtr<ExpansionMatrix> em = em_space->MakeNewExpansionMatrix();
    SmartPtr<DenseVectorSpace> expanded_space = new DenseVectorSpace(small_dim);
    SmartPtr<DenseVector> new_vec = expanded_space->MakeNewDenseVector();

    em->TransMultVector(1.0,*original,0.0,*new_vec);

    return new_vec;
  }
  printf("\nLinearizeKKT::shrink(): Unknown ERROR.");
  return NULL;
}

SmartPtr<const DenseVector> LinearizeKKT::computeDeltaP(const IntervalInfoSet& intervals,const Index& interval) const
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

SplitDecision LinearizeKKT::decideInterval(const std::vector<SplitApproximation>& approximates) const
{
  /* waiting to be implemeted*/
}

SplitApproximation LinearizeKKT::applyAlgorithmOnInterval(SmartPtr<IpoptApplication> app, const Index& interval)
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
  SmartPtr<ShiftVector> rhs, rhs_scaled, lhs;

  // start algorithm on each existing intervalID
  Index s_count = 0;
  for (Index j=0;j<intervalIDs.size();j++) {
    if (intervalIDs[j]==interval) {

      // assign private attributes depending on chosen interval
      //      if (assignIntAndParaDepAttr(interval,j)) {

    // split interval indices into to be shifted interval entries and remainder
  if (!splitIntervalIndices(x_intervalIDs_,shift_x_indices_,interval))
    printf("\nLinearizeKKT::assignIntAndParaDepAttr(): ERROR: unable to split x_intervalIDs_");
  if (!splitIntervalIndices(c_intervalIDs_,shift_c_indices_,interval))
    printf("\nLinearizeKKT::assignIntAndParaDepAttr(): ERROR: unable to split c_intervalIDs_");
  if (!splitIntervalIndices(d_intervalIDs_,shift_d_indices_,interval))
    printf("\nLinearizeKKT::assignIntAndParaDepAttr(): ERROR: unable to split d_intervalIDs_");

  // get different rhs entries for this column
  x_i_ = extractColumn(Wp_,j);
  y_c_i_ = extractColumn(Ap_,j);
  y_d_i_ = extractColumn(Bp_,j);

  // assign rhs subvectors
  u_i_ = shrink(x_i_,u_indices_);

	// create ShiftVector for computation calls
	SmartPtr<IteratesVector> rhs_top = app->IpoptDataObject()->curr()->MakeNewIteratesVector();
	SmartPtr<DenseVector> rhs_x,rhs_y_c,rhs_y_d;
	rhs_top->Set(0.0);
	rhs_top->Set_x_NonConst(*x_i_);
	rhs_top->Set_y_c_NonConst(*y_c_i_);
	rhs_top->Set_y_d_NonConst(*y_d_i_);
	rhs_x = shrink(x_i_,shift_x_indices_);
	rhs_y_c = shrink(y_c_i_,shift_c_indices_);
	rhs_y_d = shrink(y_d_i_,shift_d_indices_);
	rhs = new ShiftVector(rhs_top,rhs_x,rhs_y_c,rhs_y_d);

	// compute solution of the GMRES-Step
	const Number steplength = computeRMultS(app,rhs)*computeSMultS(app,rhs);

	if (steplength) {
	  rhs_scaled = rhs;
	  rhs_scaled->Scal(steplength);

	} else
	  printf("\nLinearizeKKT::applyAlgorithmOnInterval(): ERROR: Step-length is 0.");

	// undo conditionning step to get the approximation of the postsplit sensitivity column
	lhs = undoConditionning(app,rhs_scaled);
	z_postsplit = lhs->getDVector();

	// add the resulting sense vector to matrix ....
	sense->SetVector(s_count,*z_postsplit);
	s_count++;
    }
  }
  // get the parameter perturbation for this split
  deltap = computeDeltaP(intervals_,interval);

  // calculate absolute sensitivity values
  sense->MultVector(1.0,*deltap,0.0,*res_sense);

  // printf("\n\n");
  // res_sense->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG,"result_f");
  // printf("\n\n");

  SplitApproximation retval;

  return retval;
}

/** given an Index-type vector containing a list of intervalIDs, loads indices matching the given interval into the shift_indices part, the rest into static_indices **/
bool LinearizeKKT::splitIntervalIndices(const std::vector<Index>& intervalIDs,std::vector<Index>& shift_indices,const Index& interval)
{
  // make sure nothings in the way of push_backs
  shift_indices.clear();

  assert(intervalIDs.size());

  // cycle through intervalIDs and assign the indices
  for (int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i] == interval)
      shift_indices.push_back(i);
  }
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
  SplitAlgorithm* splitter = new LinearizeKKT(app);
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

Number Abs(Number value) // use fabs(double)
{
  Number retval =0;
  if (value<0)
    retval = (-value);
  else
    retval = value;
  return retval;
}

Number LinearizeKKT::computeRMultS(SmartPtr<IpoptApplication> app,SmartPtr<ShiftVector> target) const
{
  Number retval =0;
  SmartPtr<ShiftVector> S = computeS(app,target);
  retval+= target->Dot(*S);
  printf("\ncomputeRMultS: retval is %e", retval);
  return retval;
}

Number LinearizeKKT::computeSMultS(SmartPtr<IpoptApplication> app,SmartPtr<ShiftVector> target) const
{
  Number retval =0;

  SmartPtr<ShiftVector> S = computeS(app,target);
  retval+= S->Dot(*S);
  printf("\ncomputeSMultS: retval is %e", retval);
  return retval;
}

/* SmartPtr<const DenseVector> LinearizeKKT::computeR() const
{
// not needed for x_0 = 0 and 0 steps.
} */

SmartPtr<ShiftVector> LinearizeKKT::computeS(SmartPtr<IpoptApplication> app,SmartPtr<ShiftVector> target) const
{
  // init retval parts
  SmartPtr<DenseVector> retval_x = target->x()->MakeNewDenseVector();
  SmartPtr<DenseVector> retval_u = shrink(dynamic_cast<DenseVector*>(target->top()->x()->MakeNewCopy()),u_indices_);
  SmartPtr<DenseVector> retval_y_c = target->y_c()->MakeNewDenseVector();
  SmartPtr<DenseVector> retval_y_d = target->y_d()->MakeNewDenseVector();

  // get (W_i,shift times rhs_i) with i= (shift_h, u)
  SmartPtr<DenseVector> tmp_rhs = expand(target->x(),shift_x_indices_,W_->NCols());
  SmartPtr<DenseVector> extractor = tmp_rhs->MakeNewDenseVector();
  W_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  // map results back to get u and shift part seperately for the appropriate parts of S
  SmartPtr<DenseVector> extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,0.0);
  SmartPtr<DenseVector> extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,0.0);

  //get (A_shift^T times rhs_i) parts with i= (u, shift_c)
  tmp_rhs = expand(target->y_c(),shift_c_indices_,A_->NRows());
  SmartPtr<DenseVectorSpace> extractor_space = new DenseVectorSpace(A_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  A_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map results back to get u and shift part seperately for the appropriate parts of S
  extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,1.0);
  extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,1.0);

  //get (B_shift^T times rhs_i) parts with i= (u, shift_d)
  tmp_rhs = expand(target->y_d(),shift_d_indices_,B_->NRows());
  extractor_space = new DenseVectorSpace(B_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  B_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map result back to get u and shift part seperately for the appropriate parts of S
  extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,1.0);
  extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,1.0);

  //get (A_shift times rhs_shift_c) part
  tmp_rhs = expand(target->x(),shift_x_indices_,A_->NCols());
  extractor = tmp_rhs->MakeNewDenseVector();
  A_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map result back to the appropriate parts of S
  retval_y_c = shrink(extractor,shift_c_indices_);

  //get (B_shift times rhs_shift_d) part
  tmp_rhs = expand(target->x(),shift_x_indices_,B_->NCols());
  extractor = tmp_rhs->MakeNewDenseVector();
  B_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map result back to the appropriate parts of S
  retval_y_d = shrink(extractor,shift_d_indices_);

  // copy top values over from target, at the same time add D-multiplication to shifted part
  SmartPtr<ShiftVector> retval = computeDMultYrhs(app,target);
  //  printf("\ncomputeS() just called computeDMultYrhs()!");

  // map and insert u part back into top part
  SmartPtr<DenseVector> top_x = dynamic_cast<DenseVector*>(target->top()->x()->MakeNewCopy());
  SmartPtr<DenseVector> top_u = expand(retval_u,u_indices_,top_x->Dim());
  top_x->AddOneVector(1.0,*top_u,1.0);
  SmartPtr<IteratesVector> retval_top = target->top()->MakeNewIteratesVector();
  retval_top->Set(0.0);
  retval_top->Set_x_NonConst(*top_x);
  retval_top->Set_y_c(*target->top()->y_c());
  retval_top->Set_y_d(*target->top()->y_d());

  // insert calculated parts into retval
  SmartPtr<ShiftVector> retpart = new ShiftVector(retval_top,retval_x,retval_y_c,retval_y_d);
  retval->AddOneVector(1.0,*retpart,1.0);

  return retval;
}

SmartPtr<ShiftVector> LinearizeKKT::computeDMultYrhs(SmartPtr<IpoptApplication> app,SmartPtr<ShiftVector> target) const
{
  // get KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

  // do actual backsolve
  SmartPtr<IteratesVector> preretval = target->top()->MakeNewIteratesVector();
  pd_solver->Solve(1.0, 0.0, *target->top(), *preretval);

  // get (K^-1 times rhs) u part and eliminate non-u part elements
  SmartPtr<DenseVector> retval_x = dynamic_cast<DenseVector*>(preretval->x()->MakeNewCopy());
  retval_x = shrink(retval_x,u_indices_);
  retval_x = expand(retval_x,u_indices_,W_->NCols());
  SmartPtr<DenseVector> retval_c;  // = preretval->y_c();
  SmartPtr<DenseVector> retval_d;  // = preretval->y_d();

  // multiply with appropriate lhs parts
  SmartPtr<DenseVectorSpace> retval_nc_space = new DenseVectorSpace(W_->NRows());
  SmartPtr<DenseVector> retval_x_nonconst = retval_nc_space->MakeNewDenseVector();
  W_->MultVector(1.0,*retval_x,0.0,*retval_x_nonconst);
  retval_nc_space = new DenseVectorSpace(A_->NRows());
  SmartPtr<DenseVector> retval_c_nonconst = retval_nc_space->MakeNewDenseVector();
  A_->MultVector(1.0,*retval_x,0.0,*retval_c_nonconst);
  retval_nc_space = new DenseVectorSpace(B_->NRows());
  SmartPtr<DenseVector> retval_d_nonconst = retval_nc_space->MakeNewDenseVector();
  B_->MultVector(1.0,*retval_x,0.0,*retval_d_nonconst);

  // reduce to the desired quantities
  retval_x = shrink(retval_x_nonconst,shift_x_indices_);
  retval_c = shrink(retval_c_nonconst,shift_c_indices_);
  retval_d = shrink(retval_d_nonconst,shift_d_indices_);

  //setup return value
  SmartPtr<ShiftVector> retval = new ShiftVector(target->top()->MakeNewIteratesVectorCopy(),retval_x,retval_c,retval_d);

  return retval;
}

SmartPtr<ShiftVector> LinearizeKKT::undoConditionning(SmartPtr<IpoptApplication> app, SmartPtr<ShiftVector> target) const
{
  SmartPtr<IteratesVector> retval_top = target->top()->MakeNewIteratesVectorCopy();


  // get KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

  // do actual backsolve
  pd_solver->Solve(1.0, 0.0, *target->top(), *retval_top);
  SmartPtr<const Vector> retval_x = retval_top->x()->MakeNewCopy();
  SmartPtr<const Vector> retval_c = retval_top->y_c()->MakeNewCopy();
  SmartPtr<const Vector> retval_d = retval_top->y_d()->MakeNewCopy();

  retval_top->Set(0.0);
  retval_top->Set_x(*retval_x);
  retval_top->Set_y_c(*retval_c);
  retval_top->Set_y_d(*retval_d);

  // set up retval
  SmartPtr<ShiftVector> retval = new ShiftVector(retval_top,target->x(),target->y_c(),target->y_d());


  return retval;
}
