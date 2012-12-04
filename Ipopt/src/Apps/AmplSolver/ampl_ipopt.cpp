// Coyright (C) 2004, 2009 International Business Machines and others.
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

// for fabs and sqrt
# include <math.h>

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
class IntervalInfo
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
  Index size() const;

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
/*
class ShiftSolver : public ReferencedObject
{
public:
  virtual SmartPtr<DenseVector> solveShiftedSystem(const ShiftVector& x0,const ShiftVector& rhs) const = 0;
};

class GMRES : public ShiftSolver
{
SmartPtr<DenseVector> solveShiftedSystem(const ShiftVector& x0,const ShiftVector& rhs) const;
};

class MINRES : public ShiftSolver
{
SmartPtr<DenseVector> solveShiftedSystem(const ShiftVector& x0,const ShiftVector& rhs) const;
};

ShiftSolver* assignShiftSolver(SmartPtr<OptionsList> options);
*/
class BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const bool& force_obi=0) const = 0;
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const std::vector<Index>& indices,
							   const bool& force_obi=0) const = 0;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const = 0;
};

BranchingCriterion* assignBranchingMethod(SmartPtr<OptionsList> options);

class RandomBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const bool& force_obi=0) const;
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const std::vector<Index>& indices,
							   const bool& force_obi=0) const;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const;
};

class LargerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const bool& force_obi=0) const;
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const std::vector<Index>& indices,
							   const bool& force_obi=0) const;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const;
};

class SmallerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const bool& force_obi=0) const;
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const std::vector<Index>& indices,
							   const bool& force_obi=0) const;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const;
};

class AbsoluteLargerBranching : public BranchingCriterion
{
public:
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const bool& force_obi=0) const;
  virtual std::vector<SplitChoice> branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,
							   SmartPtr<OptionsList> options,
							   const IntervalInfoSet& intervals,
							   const std::vector<Index>& indices,
							   const bool& force_obi=0) const;
  virtual SplitChoice branchInterval(const std::vector<SplitChoice>& splitchoices) const;
};

class Intervaluation
{
public:
  virtual SplitChoice intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals) const= 0;
};
Intervaluation* assignIntervaluationMethod(SmartPtr<OptionsList> options);

class OneBoundIntervaluation : public Intervaluation
{
public:
  SplitChoice intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals) const;
};

class BothBoundIntervaluation : public Intervaluation
{
public:
  SplitChoice intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals) const;
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
  SmartPtr<DenseVector> applyGMRESOnInterval(SmartPtr<IpoptApplication> app,const Index& interval, const Index& parameter);

  /*handling and manipulating data specifically*/
  bool splitIntervalIndices(const std::vector<Index>& intervalIDs,std::vector<Index>& shift_indices,const Index& interval);
  SmartPtr<DenseVector> extractColumn(SmartPtr<const Matrix> original,const Index& column) const;
  SmartPtr<DenseVector> expand(SmartPtr<DenseVector> original,const std::vector<Index>& indices, const Index& large_dim) const;
  SmartPtr<DenseVector> expand(SmartPtr<DenseVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<DenseVector> expand(SmartPtr<IteratesVector> original, const Index& large_dim, const Index& start_idx) const;
  SmartPtr<DenseVector> shrink(SmartPtr<DenseVector> original,const std::vector<Index>& indices) const;
  SmartPtr<MultiVectorMatrix> computeVMultEye(SmartPtr<const DenseVector> target) const;
  /*algorithm specific components*/
  void assignGMRESOptions(SmartPtr<OptionsList> options,Number& tol, Index& n_max, Index& n_rest) const;
  SmartPtr<ShiftVector> computeGMRES(SmartPtr<ShiftVector>b,SmartPtr<ShiftVector>x0,const Number& tol =1.0e-6, const Index& n_max=5, const Index& n_rest=0);
  SmartPtr<ShiftVector> computeAMultVector(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeBMultVector(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeCMultVector(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeDMultVector(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeAMVminres(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeBMVminres(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeCMVminres(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeDMVminres(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeKSymMultVector(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeKIMultVector(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeKPIMultVector(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computePMultVector(SmartPtr<ShiftVector> target) const;
  SmartPtr<ShiftVector> computeR(SmartPtr<ShiftVector> b, SmartPtr<ShiftVector>x) const;
  SmartPtr<ShiftVector> computeV(SmartPtr<ShiftVector> r) const;
  SmartPtr<const DenseVector> computeDeltaP(const IntervalInfoSet& intervals,const Index& interval, const Index& parameterID) const;
  SmartPtr<DenseVector> applyMINRESOnInterval(SmartPtr<IpoptApplication> app,const Index& interval, const Index& parameter);
  SmartPtr<ShiftVector> computeMINRES(SmartPtr<ShiftVector>b,SmartPtr<ShiftVector>x0,const Number& tol =1.0e-6, const Index& n_max=5, const Index& n_rest=0) const;
  void convertBackToUnsymmetric(SmartPtr<ShiftVector> target) const;
  std::vector<Index> collectIntMetaData() const;
  /*linear system solution algorithm test functions*/
  void testMINRES();
  void testGMRES();

private:
  //  SmartPtr<ShiftSolver> solver_;
  //original Ipopt data

  SmartPtr<const DenseVector> x_;
  SmartPtr<const DenseVector> s_;
  SmartPtr<const DenseVector> y_c_;
  SmartPtr<const DenseVector> y_d_;
  SmartPtr<const DenseVector> z_L_;
  SmartPtr<const DenseVector> z_U_;
  SmartPtr<const DenseVector> v_L_;
  SmartPtr<const DenseVector> v_U_;

  SmartPtr<const Matrix> Wp_;
  SmartPtr<const Matrix> Ap_;
  SmartPtr<const Matrix> Bp_;
  SmartPtr<const SymMatrix> W_;
  SmartPtr<const Matrix> A_;
  SmartPtr<const Matrix> B_;
  SmartPtr<PDSystemSolver> KKT_;
  SmartPtr<Journalist> jnl_;

  SmartPtr<const Vector> x_L_;    // Lower bounds on x
  SmartPtr<const Matrix> Px_L_;   // Permutation matrix (x_L_ -> x)
  SmartPtr<const Vector> x_U_;    // Upper bounds on x
  SmartPtr<const Matrix> Px_U_;   // Permutation matrix (x_U_ -> x)
  SmartPtr<const Vector> d_L_;    // Lower bounds on d
  SmartPtr<const Matrix> Pd_L_;   // Permutation matrix (d_L_ -> d)
  SmartPtr<const Vector> d_U_;    // Upper bounds on d
  SmartPtr<const Matrix> Pd_U_;   // Permutation matrix (d_U_ -> d)

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
  Index n_ul_;
  Index n_uu_;

  SmartPtr<DenseVector> z_L_DIV_Sl_;
  SmartPtr<DenseVector> z_U_DIV_Sl_;
  SmartPtr<DenseVector> v_L_DIV_Sl_;
  SmartPtr<DenseVector> v_U_DIV_Sl_;
  SmartPtr<DenseVector> Sl_x_L_;
  SmartPtr<DenseVector> Sl_x_U_;
  SmartPtr<DenseVector> Sl_s_L_;
  SmartPtr<DenseVector> Sl_s_U_;

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

class SplitIntAtBound : public SplitAlgorithm
{
  SplitDecision applySplitAlgorithm(SmartPtr<IpoptApplication> app);
};


/*forward declaration*/
//class ShiftVectorSpace;

class ShiftVector : public ReferencedObject
{
public:
  ShiftVector(SmartPtr<IteratesVector> top,
	      SmartPtr<DenseVector> x,
	      SmartPtr<DenseVector> s,
	      SmartPtr<DenseVector>y_c,
	      SmartPtr<DenseVector> y_d,
	      SmartPtr<DenseVector>z_L,
	      SmartPtr<DenseVector>z_U,
	      SmartPtr<DenseVector>v_L,
	      SmartPtr<DenseVector>v_U);
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
  //  Number Dot(const Vector &x) const;
  Number Dot(const ShiftVector &x) const;
  Index Dim() const;
  void Set(const Number& alpha);
  Number Nrm2() const;
  void AddOneVector(const Number& a,const ShiftVector& v1,const Number& c);
  void AddTwoVectors(const Number& a,const ShiftVector& v1,
		     const Number& b,const ShiftVector& v2,const Number& c);
  SmartPtr<ShiftVector> MakeNewShiftVector() const;
  //  SmartPtr<Vector> MakeNew() const; // disabled: no Vector child yet
  SmartPtr<IteratesVector>top() const;
  SmartPtr<DenseVector>x() const;
  SmartPtr<DenseVector>s() const;
  SmartPtr<DenseVector>y_c() const;
  SmartPtr<DenseVector>y_d() const;
  SmartPtr<DenseVector>z_L() const;
  SmartPtr<DenseVector>z_U() const;
  SmartPtr<DenseVector>v_L() const;
  SmartPtr<DenseVector>v_U() const;
  void appendValues(SmartPtr<const DenseVector>x,std::vector<Number>& values) const;
  void Set_x_top(const DenseVector& x);
  void Set_s_top(const DenseVector& s);
  void Set_y_c_top(const DenseVector& y_c);
  void Set_y_d_top(const DenseVector& y_d);
  void Set_z_L_top(const DenseVector& z_L);
  void Set_z_U_top(const DenseVector& z_U);
  void Set_v_L_top(const DenseVector& v_L);
  void Set_v_U_top(const DenseVector& v_U);
  void Set_x(const DenseVector& x);
  void Set_s(const DenseVector& s);
  void Set_y_c(const DenseVector& y_c);
  void Set_y_d(const DenseVector& y_d);
  void Set_z_L(const DenseVector& z_L);
  void Set_z_U(const DenseVector& z_U);
  void Set_v_L(const DenseVector& v_L);
  void Set_v_U(const DenseVector& v_U);

private:
  ShiftVector();
  SmartPtr<DenseVectorSpace> x_space;
  SmartPtr<DenseVectorSpace> y_c_space;
  SmartPtr<DenseVectorSpace> y_d_space;
  SmartPtr<const IteratesVectorSpace> t_space;

  SmartPtr<IteratesVector> top_;
  SmartPtr<DenseVector> x_;
  SmartPtr<DenseVector> s_;
  SmartPtr<DenseVector> y_c_;
  SmartPtr<DenseVector> y_d_;
  SmartPtr<DenseVector> z_L_;
  SmartPtr<DenseVector> z_U_;
  SmartPtr<DenseVector> v_L_;
  SmartPtr<DenseVector> v_U_;
  //  SmartPtr<ShiftVectorSpace> owner_space_;
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
  unsigned int tmp_index;
  unsigned int i=0;

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
  if ( (unsigned) intindex<indexvec_.size())
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
  if ( (unsigned) intindex<intervalIDvec_.size())
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
  for (unsigned int i=0;i<intervalIDvec_.size();i++)
    if (intervalIDvec_[i]<*tmp_new_ID)
      *tmp_new_ID=intervalIDvec_[i];
  tmp_IDs.push_back(*tmp_new_ID);

  //assign follow up values
  while (tmp_IDs.size()<(unsigned) nints) {
    *tmp_ID = *tmp_new_ID;
    for (unsigned int i=0;i<intervalIDvec_.size();i++) {
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
  if ((unsigned)paraindex<parameterIDvec_.size())
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
  for (unsigned int i=0;i<parameterIDvec_.size();i++)
    if (parameterIDvec_[i]<*tmp_new_ID)
      *tmp_new_ID=parameterIDvec_[i];
  tmp_IDs.push_back(*tmp_new_ID);

  //assign follow up values
  while (tmp_IDs.size()<(unsigned)npars) {
    *tmp_ID = *tmp_new_ID;
    for (unsigned int i=0;i<parameterIDvec_.size();i++) {
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
  if ((unsigned)isupperindex<is_uppervec_.size())
    return is_uppervec_[isupperindex];
  else {
    printf("\nAmplTNLP.cpp: ERROR: IntervalInfoSet::isUpper() call with out of range index!\n");
    printf("\nsize of is_uppervec_: %d \n", int(is_uppervec_.size()));
    return 0;
  }
}

Index IntervalInfoSet::getOtherBndIdx(const Index& boundindex) const
{
  for (unsigned int i=0;i<intinfovec_.size();i++){
    if (parameterIDvec_[boundindex]==parameterIDvec_[i] && intervalIDvec_[boundindex]==intervalIDvec_[i] && is_uppervec_[boundindex]!=is_uppervec_[i]){
      return indexvec_[i];
    }
  }
  return 0;
}

Index IntervalInfoSet::getParameterCount() const
{
  Index tmp_count =0;
  for (unsigned int i=0;i<parameterIDvec_.size();i++) {
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
  for (unsigned int i=0;i<intervalIDvec_.size();i++) {
    if (i==0)
      tmp_count = intervalIDvec_[i];
    if (intervalIDvec_[i]>tmp_count)
      tmp_count = intervalIDvec_[i];
  }
  return tmp_count;
}

void IntervalInfoSet::printSet() const
{
  for (unsigned int i=0; i<intinfovec_.size();i++){
    printf("\n\nIntervalInfoSet Eintrag %d:", i);
    intinfovec_[i].printSet();
  }
}

Index IntervalInfoSet::size() const
{
  return int(intinfovec_.size());
}

ShiftVector::ShiftVector(SmartPtr<IteratesVector> top, SmartPtr<DenseVector> x, SmartPtr<DenseVector> s, SmartPtr<DenseVector>y_c, SmartPtr<DenseVector> y_d, SmartPtr<DenseVector> z_L, SmartPtr<DenseVector> z_U, SmartPtr<DenseVector> v_L, SmartPtr<DenseVector> v_U)
  :
  top_(top),
  x_(x),
  s_(s),
  y_c_(y_c),
  y_d_(y_d),
  z_L_(z_L),
  z_U_(z_U),
  v_L_(v_L),
  v_U_(v_U)
{
  assert((GetRawPtr(top)) && (GetRawPtr(x)) && (GetRawPtr(s)) && (GetRawPtr(y_c)) && (GetRawPtr(y_d))
	 && (GetRawPtr(z_L)) && (GetRawPtr(z_U)) && (GetRawPtr(v_L)) && (GetRawPtr(v_U)));
}

ShiftVector::ShiftVector(const ShiftVector& rhs)
  :
  top_(rhs.top()->MakeNewIteratesVectorCopy())
{
  Vector* xvec = rhs.x()->MakeNewCopy();
  Vector* svec = rhs.s()->MakeNewCopy();
  Vector* ycvec = rhs.y_c()->MakeNewCopy();
  Vector* ydvec = rhs.y_d()->MakeNewCopy();
  Vector* zlvec = rhs.z_L()->MakeNewCopy();
  Vector* zuvec = rhs.z_U()->MakeNewCopy();
  Vector* vlvec = rhs.v_L()->MakeNewCopy();
  Vector* vuvec = rhs.v_U()->MakeNewCopy();
  x_ = (dynamic_cast<DenseVector*>(xvec));
  s_ = (dynamic_cast<DenseVector*>(svec));
  y_c_ = (dynamic_cast<DenseVector*>(ycvec));
  y_d_ = (dynamic_cast<DenseVector*>(ydvec));
  z_L_ = (dynamic_cast<DenseVector*>(zlvec));
  z_U_ = (dynamic_cast<DenseVector*>(zuvec));
  v_L_ = (dynamic_cast<DenseVector*>(vlvec));
  v_U_ = (dynamic_cast<DenseVector*>(vuvec));
}

ShiftVector& ShiftVector::operator=(const ShiftVector &rhs)
{
  if (this!= &rhs)
    {
      top_ = rhs.top()->MakeNewIteratesVectorCopy();
      Vector* xvec = rhs.x()->MakeNewCopy();
      x_ = dynamic_cast<DenseVector*>(xvec);
      Vector* svec = rhs.s()->MakeNewCopy();
      s_ = dynamic_cast<DenseVector*>(svec);
      Vector* ycvec = rhs.y_c()->MakeNewCopy();
      y_c_ = dynamic_cast<DenseVector*>(ycvec);
      Vector* ydvec = rhs.y_d()->MakeNewCopy();
      y_d_ = dynamic_cast<DenseVector*>(ydvec);
      Vector* zlvec = rhs.z_L()->MakeNewCopy();
      z_L_ = dynamic_cast<DenseVector*>(zlvec);
      Vector* zuvec = rhs.z_U()->MakeNewCopy();
      z_U_ = dynamic_cast<DenseVector*>(zuvec);
      Vector* vlvec = rhs.v_L()->MakeNewCopy();
      v_L_ = dynamic_cast<DenseVector*>(vlvec);
      Vector* vuvec = rhs.v_U()->MakeNewCopy();
      v_U_ = dynamic_cast<DenseVector*>(vuvec);
    }
  return *this;
}

void ShiftVector::appendValues(SmartPtr<const DenseVector>x,std::vector<Number>& values) const
{
  assert (GetRawPtr(x));

  const Index dim = x->Dim();

  if (x->IsHomogeneous()) {
    const Number val = x->Scalar();
    for (int i=0;i<dim;i++) {
      values.push_back(val);
    }
  } else {
    const Number* vals = x->Values();
    for (int i=0;i<dim;i++) {
      values.push_back(vals[i]);
    }
  }
}

SmartPtr<DenseVector> ShiftVector::getDVector() const
{
  SmartPtr<DenseVectorSpace> retval_space = new DenseVectorSpace(Dim());
  SmartPtr<DenseVector> retval = retval_space->MakeNewDenseVector();

  Number* exp_values = new Number[Dim()];
  std::vector<Number> rvalues;

  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(top_->x())),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(top_->s())),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(top_->y_c())),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(top_->y_d())),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(top_->z_L())),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(top_->z_U())),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(top_->v_L())),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(top_->v_U())),rvalues);

  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(x_)),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(s_)),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(y_c_)),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(y_d_)),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(z_L_)),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(z_U_)),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(v_L_)),rvalues);
  appendValues(dynamic_cast<const DenseVector*>(GetRawPtr(v_U_)),rvalues);

  for (int i=0;i<Dim();i++) {
    exp_values[i] = rvalues[i];
  }

  retval->SetValues(exp_values);
  return retval;
}

void ShiftVector::Scal(const Number& factor)
{
  top_->Scal(factor);
  x_->Scal(factor);
  s_->Scal(factor);
  y_c_->Scal(factor);
  y_d_->Scal(factor);
  z_L_->Scal(factor);
  z_U_->Scal(factor);
  v_L_->Scal(factor);
  v_U_->Scal(factor);
}
/*
SmartPtr<ShiftVectorSpace> ShiftVector::OwnerSpace() const
{
  return owner_space_;
  } */

void ShiftVector::Print(SmartPtr<const Journalist> jnlst, EJournalLevel level, EJournalCategory category, const std::string& name, Index indent, const std::string& prefix) const
{
  printf("\ntop\n");
  top_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nx\n");
  x_->Print(jnlst,level,category,name,indent,prefix);
  printf("\ns\n");
  s_->Print(jnlst,level,category,name,indent,prefix);
  printf("\ny_c\n");
  y_c_->Print(jnlst,level,category,name,indent,prefix);
  printf("\ny_d\n");
  y_d_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nz_L\n");
  z_L_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nz_U\n");
  z_U_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nv_L\n");
  v_L_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nv_U\n");
  v_U_->Print(jnlst,level,category,name,indent,prefix);
}

void ShiftVector::Print(const Journalist& jnlst, EJournalLevel level, EJournalCategory category, const std::string& name, Index indent, const std::string& prefix) const
{
  printf("\ntop\n");
  top_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nx\n");
  x_->Print(jnlst,level,category,name,indent,prefix);
  printf("\ns\n");
  s_->Print(jnlst,level,category,name,indent,prefix);
  printf("\ny_c\n");
  y_c_->Print(jnlst,level,category,name,indent,prefix);
  printf("\ny_d\n");
  y_d_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nz_L\n");
  z_L_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nz_U\n");
  z_U_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nv_L\n");
  v_L_->Print(jnlst,level,category,name,indent,prefix);
  printf("\nv_U\n");
  v_U_->Print(jnlst,level,category,name,indent,prefix);
}

Number ShiftVector::Dot(const ShiftVector &x) const
{
  assert(x.Dim()==Dim());
  Number retval = 0.0;
  retval += x.top()->Dot(*top_);
  retval += x.x()->Dot(*x_);
  retval += x.s()->Dot(*s_);
  retval += x.y_c()->Dot(*y_c_);
  retval += x.y_d()->Dot(*y_d_);
  retval += x.z_L()->Dot(*z_L_);
  retval += x.z_U()->Dot(*z_U_);
  retval += x.v_L()->Dot(*v_L_);
  retval += x.v_U()->Dot(*v_U_);
  return retval;
}

void ShiftVector::AddOneVector(const Number& a,const ShiftVector& v1,const Number& c)
{
  AddTwoVectors(a,v1,0.0,*this,c);
}

void ShiftVector::AddTwoVectors(const Number& a,const ShiftVector& v1,const Number& b,const ShiftVector& v2,const Number& c)
{
  top_->AddTwoVectors(a,*v1.top(),b,*v2.top(),c);
  x_->AddTwoVectors(a,*v1.x(),b,*v2.x(),c);
  s_->AddTwoVectors(a,*v1.s(),b,*v2.s(),c);
  y_c_->AddTwoVectors(a,*v1.y_c(),b,*v2.y_c(),c);
  y_d_->AddTwoVectors(a,*v1.y_d(),b,*v2.y_d(),c);
  z_L_->AddTwoVectors(a,*v1.z_L(),b,*v2.z_L(),c);
  z_U_->AddTwoVectors(a,*v1.z_U(),b,*v2.z_U(),c);
  v_L_->AddTwoVectors(a,*v1.v_L(),b,*v2.v_L(),c);
  v_U_->AddTwoVectors(a,*v1.v_U(),b,*v2.v_U(),c);
}

Index ShiftVector::Dim() const
{
  Index retval = 0;
  retval += top_->Dim();
  retval += x_->Dim();
  retval += s_->Dim();
  retval += y_c_->Dim();
  retval += y_d_->Dim();
  retval += z_L_->Dim();
  retval += z_U_->Dim();
  retval += v_L_->Dim();
  retval += v_U_->Dim();
  return retval;
}

void ShiftVector::Set(const Number& alpha)
{
  dynamic_cast<DenseVector*>(GetRawPtr(top_->x_NonConst()))->Set(alpha);
  dynamic_cast<DenseVector*>(GetRawPtr(top_->s_NonConst()))->Set(alpha);
  dynamic_cast<DenseVector*>(GetRawPtr(top_->y_c_NonConst()))->Set(alpha);
  dynamic_cast<DenseVector*>(GetRawPtr(top_->y_d_NonConst()))->Set(alpha);
  dynamic_cast<DenseVector*>(GetRawPtr(top_->z_L_NonConst()))->Set(alpha);
  dynamic_cast<DenseVector*>(GetRawPtr(top_->z_U_NonConst()))->Set(alpha);
  dynamic_cast<DenseVector*>(GetRawPtr(top_->v_L_NonConst()))->Set(alpha);
  dynamic_cast<DenseVector*>(GetRawPtr(top_->v_U_NonConst()))->Set(alpha);
  x_->Set(alpha);
  s_->Set(alpha);
  y_c_->Set(alpha);
  y_d_->Set(alpha);
  z_L_->Set(alpha);
  z_U_->Set(alpha);
  v_L_->Set(alpha);
  v_U_->Set(alpha);
}

Number ShiftVector::Nrm2() const
{
  Number retval = top_->Dot(*top_) + x_->Dot(*x_) + y_c_->Dot(*y_c_) + y_d_->Dot(*y_d_);
  retval+=  s_->Dot(*s_) + z_L_->Dot(*z_L_) + z_U_->Dot(*z_U_) + v_L_->Dot(*v_L_) + v_U_->Dot(*v_U_);
  return sqrt(retval);
}

SmartPtr<ShiftVector> ShiftVector::MakeNewShiftVector() const
{
  SmartPtr<IteratesVector> top = dynamic_cast<IteratesVector*>(top_->MakeNew());
  SmartPtr<DenseVector> x = x_->MakeNewDenseVector();
  SmartPtr<DenseVector> s = s_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_c = y_c_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_d = y_d_->MakeNewDenseVector();
  SmartPtr<DenseVector> z_L = z_L_->MakeNewDenseVector();
  SmartPtr<DenseVector> z_U = z_U_->MakeNewDenseVector();
  SmartPtr<DenseVector> v_L = v_L_->MakeNewDenseVector();
  SmartPtr<DenseVector> v_U = v_U->MakeNewDenseVector();
  return  new ShiftVector(top,x,s,y_c,y_d,z_L,z_U,v_L,v_U);
}
/*  // disabled: no Vector child yet
SmartPtr<Vector> ShiftVector::MakeNew() const
{
  return exp_->MakeNew();
}
*/

SmartPtr<IteratesVector> ShiftVector::top() const
{
  return top_;
}

SmartPtr<DenseVector> ShiftVector::x() const
{
  return x_;
}

SmartPtr<DenseVector> ShiftVector::s() const
{
  return s_;
}

SmartPtr<DenseVector> ShiftVector::y_c() const
{
  return y_c_;
}

SmartPtr<DenseVector> ShiftVector::y_d() const
{
  return y_d_;
}

SmartPtr<DenseVector> ShiftVector::z_L() const
{
  return z_L_;
}

SmartPtr<DenseVector> ShiftVector::z_U() const
{
  return z_U_;
}

SmartPtr<DenseVector> ShiftVector::v_L() const
{
  return v_L_;
}

SmartPtr<DenseVector> ShiftVector::v_U() const
{
  return v_U_;
}

void ShiftVector::Set_x_top(const DenseVector& x)
{
  top_->Set_x(x);
}

void ShiftVector::Set_s_top(const DenseVector& s)
{
  top_->Set_s(s);
}

void ShiftVector::Set_y_c_top(const DenseVector& y_c)
{
  top_->Set_y_c(y_c);
}

void ShiftVector::Set_y_d_top(const DenseVector& y_d)
{
  top_->Set_y_d(y_d);
}

void ShiftVector::Set_z_L_top(const DenseVector& z_L)
{
  top_->Set_z_L(z_L);
}

void ShiftVector::Set_z_U_top(const DenseVector& z_U)
{
  top_->Set_z_U(z_U);
}

void ShiftVector::Set_v_L_top(const DenseVector& v_L)
{
  top_->Set_v_L(v_L);
}

void ShiftVector::Set_v_U_top(const DenseVector& v_U)
{
  top_->Set_v_U(v_U);
}

void ShiftVector::Set_x(const DenseVector& x)
{
  if (x.IsHomogeneous())
    x_->Set(x.Scalar());
  else {
    const Number* vals = x.Values();
    x_->SetValues(vals);
  }
}

void ShiftVector::Set_s(const DenseVector& s)
{
  if (s.IsHomogeneous())
    s_->Set(s.Scalar());
  else {
    const Number* vals = s.Values();
    s_->SetValues(vals);
  }
}

void ShiftVector::Set_y_c(const DenseVector& y_c)
{
  if (y_c.IsHomogeneous())
    y_c_->Set(y_c.Scalar());
  else {
    const Number* vals = y_c.Values();
    y_c_->SetValues(vals);
  }
}

void ShiftVector::Set_y_d(const DenseVector& y_d)
{
  if (y_d.IsHomogeneous())
    y_d_->Set(y_d.Scalar());
  else {
    const Number* vals = y_d.Values();
    y_d_->SetValues(vals);
  }
}

void ShiftVector::Set_z_L(const DenseVector& z_L)
{
  if (z_L.IsHomogeneous())
    z_L_->Set(z_L.Scalar());
  else {
    const Number* vals = z_L.Values();
    z_L_->SetValues(vals);
  }
}

void ShiftVector::Set_z_U(const DenseVector& z_U)
{
  if (z_U.IsHomogeneous())
    z_U_->Set(z_U.Scalar());
  else {
    const Number* vals = z_U.Values();
    z_U_->SetValues(vals);
  }
}

void ShiftVector::Set_v_L(const DenseVector& v_L)
{
  if (v_L.IsHomogeneous())
    v_L_->Set(v_L.Scalar());
  else {
    const Number* vals = v_L.Values();
    v_L_->SetValues(vals);
  }
}

void ShiftVector::Set_v_U(const DenseVector& v_U)
{
  if (v_U.IsHomogeneous())
    v_U_->Set(v_U.Scalar());
  else {
    const Number* vals = v_U.Values();
    v_U_->SetValues(vals);
  }
}

/*  else {
    const Number* t_y_c_values = dynamic_cast<const DenseVector*>(GetRawPtr(top->y_c()))->Values();
    if (!t_y_c_values)
      printf("\nShiftVector::ShiftVector(): Error: t_y_c_values is NULL.");
    else {
      Number* t_y_c__values = new Number[t_y_c_->Dim()];
      std::copy(t_y_c_values, t_y_c_values+t_y_c_->Dim(),t_y_c__values);
      t_y_c_->SetValues(t_y_c__values);
    }
  }
  t_y_d_space = new DenseVectorSpace(top->y_d()->Dim());
  SmartPtr<DenseVector> t_y_d_ = t_y_d_space->MakeNewDenseVector();

if (dynamic_cast<const DenseVector*>(GetRawPtr(top->y_d()))->IsHomogeneous())
    t_y_d_->Set(dynamic_cast<const DenseVector*>(GetRawPtr(top->y_d()))->Scalar());
  else {
    const Number* t_y_d_values = dynamic_cast<const DenseVector*>(GetRawPtr(top->y_d()))->Values();
    if (!t_y_d_values)
      printf("\nShiftVector::ShiftVector(): Error: t_y_d_values is NULL.");
    else {
      Number* t_y_d__values = new Number[t_y_d_->Dim()];
      std::copy(t_y_d_values, t_y_d_values+t_y_d_->Dim(),t_y_d__values);
      t_y_d_->SetValues(t_y_d__values);
    }
  }
    top_->Set_x_NonConst(*t_x_);
    top_->Set_y_c_NonConst(*t_y_c_);
    top_->Set_y_d_NonConst(*t_y_d_);
}

void ShiftVector::x(SmartPtr<DenseVector>x)
{
  x_space = new DenseVectorSpace(x->Dim());
  x_ = x_space->MakeNewDenseVector();

  if (x->IsHomogeneous())
    x_->Set(x->Scalar());
  else {
    const Number* x_values = x->Values();
    if (!x_values)
      printf("\nShiftVector::ShiftVector(): Error: x_values is NULL.");
    Number* x__values = new Number[x->Dim()];
    std::copy(x_values, x_values+x->Dim(),x__values);
    x_->SetValues(x__values);
  }
  x_dim_ = x_->Dim();
}

void ShiftVector::y_c(SmartPtr<DenseVector>y_c)
{
  y_c_space = new DenseVectorSpace(y_c->Dim());
  y_c_ = y_c_space->MakeNewDenseVector();

  if (y_c->IsHomogeneous())
    y_c_->Set(y_c->Scalar());
  else {
    const Number* y_c_values = y_c->Values();
    if (!y_c_values)
      printf("\nShiftVector::ShiftVector(): Error: y_c_values is NULL.");
    Number* y_c__values = new Number[y_c->Dim()];
    std::copy(y_c_values, y_c_values+y_c->Dim(),y_c__values);
    y_c_->SetValues(y_c__values);
  }
  y_c_dim_ = y_c_->Dim();
}

void ShiftVector::y_d(SmartPtr<DenseVector>y_d)
{
  y_d_space = new DenseVectorSpace(y_d->Dim());
  y_d_ = y_d_space->MakeNewDenseVector();

  if (y_d->IsHomogeneous())
    y_d_->Set(y_d->Scalar());
  else {
    const Number* y_d_values = y_d->Values();
    if (!y_d_values)
      printf("\nShiftVector::ShiftVector(): Error: y_d_values is NULL.");
    Number* y_d__values = new Number[y_d->Dim()];
    std::copy(y_d_values, y_d_values+y_d->Dim(),y_d__values);
    y_d_->SetValues(y_d__values);
  }
  y_d_dim_ = y_d_->Dim();
}

void ShiftVector::ownerspace(SmartPtr<ShiftVectorSpace>ospace)
{
  owner_space_ = ospace;
  dim_ = ospace->Dim();
  }

void ShiftVector::exp(SmartPtr<DenseVector> exp)
{
  exp_ = exp;
  exp_space_ = new DenseVectorSpace(exp_->Dim());
} */
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

  // disabled: no Vector child yet
SmartPtr<Vector>ShiftVectorSpace::MakeNew() const
{
  SmartPtr<IteratesVector> top = top_space_->MakeNewIteratesVector();
  SmartPtr<DenseVector> x = x_space_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_c = y_c_space_->MakeNewDenseVector();
  SmartPtr<DenseVector> y_d = y_d_space_->MakeNewDenseVector();

  return dynamic_cast<Vector*>(new ShiftVector(top,x,y_c,y_d));
}

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
    for (int j=0;j<intervals.size();j++) {
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
    for (int j=0;j<intervals.size();j++) {
      if (parameterflags[j] == parameterIDs[i] && intervalflags[j] == intervalID) {
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

std::vector<SplitChoice> RandomBranching::branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals, const std::vector<Index>& indices,const bool& force_obi) const
{
  // RandomBranching not implemented.
}

std::vector<SplitChoice> RandomBranching::branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals, const bool& force_obi) const
{
  // RandomBranching not implemented.
}


SplitChoice RandomBranching::branchInterval(const std::vector<SplitChoice>& splitchoices) const
{
  // RandomBranching not implemented.
}

std::vector<SplitChoice> LargerBranching::branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals, const std::vector<Index>& indices,const bool& force_obi) const
{
  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs = intervals.getIntervalIDs();
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater;
  if (!force_obi)
    intervaluater = assignIntervaluationMethod(options);
  else
    intervaluater = new OneBoundIntervaluation();

  assert(intervalIDs.size());
  for (unsigned int i=0;i<indices.size();i++) {
    //with only one control and only one interval, the only decision left is which parameter to split. thats done by intervaluation
    retval.push_back(intervaluater->intervaluateInterval(indices[i],intervalIDs[0],mv_sens,options,intervals));
    if (intervalIDs.size()>1) {
      //with more than one interval, the branchingalgorithm is to pick the critical interval and return its data
      for (unsigned int j=1;j<intervalIDs.size();j++) {
	tmp_split_choice = intervaluater->intervaluateInterval(indices[i],intervalIDs[j],mv_sens,options,intervals);
	if (tmp_split_choice.reason_value>retval[i].reason_value)
	  retval[i] = tmp_split_choice;
      }
    }
  }
  return retval;
}

std::vector<SplitChoice>  LargerBranching::branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals, const bool& force_obi) const
{
  // cycle through var space interval flags to identify and save control indexes
  const Index nrows = mv_sens->NRows();
  std::vector<Index> ctrl_rows;
  assert(dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->HasIntegerMetaData("intervalID"));
  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");
  for (int i=0;i<nrows;i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
  }
  return branchSensitivityMatrix(mv_sens,options,intervals,ctrl_rows,force_obi);
}

SplitChoice LargerBranching::branchInterval(const std::vector<SplitChoice>& splitchoices) const
{
  assert(splitchoices.size());
  SplitChoice retval=splitchoices[0];
  // otherwise pick the interval critical to branchmode
  for (unsigned int i=1;i<splitchoices.size();i++) {
    if (splitchoices[i].reason_value>retval.reason_value)
      retval=splitchoices[i];
  }
  return retval;
}

std::vector<SplitChoice> SmallerBranching::branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals, const std::vector<Index>& indices,const bool& force_obi) const
{
  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs = intervals.getIntervalIDs();
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater;
  if (!force_obi)
    intervaluater = assignIntervaluationMethod(options);
  else
    intervaluater = new OneBoundIntervaluation();
  assert(intervalIDs.size());
  for (unsigned int i=0;i<indices.size();i++) {
    //with only one control and only one interval, the only decision left is which parameter to split. thats done by intervaluation
    retval.push_back(intervaluater->intervaluateInterval(indices[i],intervalIDs[0],mv_sens,options,intervals));
    if (intervalIDs.size()>1) {
      //with more than one interval, the branchingalgorithm is to pick the critical interval and return its data
      for (unsigned int j=1;j<intervalIDs.size();j++) {
	tmp_split_choice = intervaluater->intervaluateInterval(indices[i],intervalIDs[j],mv_sens,options,intervals);
	if (tmp_split_choice.reason_value<retval[i].reason_value)
	  retval[i] = tmp_split_choice;
      }
    }
  }
  return retval;
}

std::vector<SplitChoice> SmallerBranching::branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals, const bool& force_obi) const
{
  // cycle through var space interval flags to identify and save control indexes
  const Index nrows = mv_sens->NRows();
  std::vector<Index> ctrl_rows;
  assert(dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->HasIntegerMetaData("intervalID"));
  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");
  for (int i=0;i<nrows;i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
  }
  return branchSensitivityMatrix(mv_sens,options,intervals,ctrl_rows,force_obi);
}

SplitChoice SmallerBranching::branchInterval(const std::vector<SplitChoice>& splitchoices) const
{
  assert(splitchoices.size());
  SplitChoice retval=splitchoices[0];
  // otherwise pick the interval critical to branchmode
  for (unsigned int i=1;i<splitchoices.size();i++) {
    if (splitchoices[i].reason_value<retval.reason_value) {
      retval=splitchoices[i];
    }
  }
  return retval;
}

std::vector<SplitChoice> AbsoluteLargerBranching::branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals, const std::vector<Index>& indices,const bool& force_obi) const
{
  std::vector<SplitChoice> retval;
  std::vector<Index> intervalIDs = intervals.getIntervalIDs();
  // cycle through controls and intervals to setup intervaluation for each control at each interval
  SplitChoice tmp_split_choice;
  Intervaluation* intervaluater;
  if (!force_obi)
    intervaluater = assignIntervaluationMethod(options);
  else
    intervaluater = new OneBoundIntervaluation();
  assert(intervalIDs.size());
  for (unsigned int i=0;i<indices.size();i++) {
    //with only one control and only one interval, the only decision left is which parameter to split. thats done by intervaluation
    retval.push_back(intervaluater->intervaluateInterval(indices[i],intervalIDs[0],mv_sens,options,intervals));
    if (intervalIDs.size()>1) {
      //with more than one interval, the branchingalgorithm is to pick the critical interval and return its data
      for (unsigned int j=1;j<intervalIDs.size();j++) {
	tmp_split_choice = intervaluater->intervaluateInterval(indices[i],intervalIDs[j],mv_sens,options,intervals);
	if (fabs(tmp_split_choice.reason_value)>fabs(retval[i].reason_value))
	  retval[i] = tmp_split_choice;
      }
    }
  }
  return retval;
}

std::vector<SplitChoice> AbsoluteLargerBranching::branchSensitivityMatrix(SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals, const bool& force_obi) const
{
  // cycle through var space interval flags to identify and save control indexes
  const Index nrows = mv_sens->NRows();
  std::vector<Index> ctrl_rows;
  assert(dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->HasIntegerMetaData("intervalID"));
  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");
  for (int i=0;i<nrows;i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
  }
  return branchSensitivityMatrix(mv_sens,options,intervals,ctrl_rows,force_obi);
}

SplitChoice AbsoluteLargerBranching::branchInterval(const std::vector<SplitChoice>& splitchoices) const
{
  assert(splitchoices.size());
  SplitChoice retval=splitchoices[0];
  // otherwise pick the interval critical to branchmode
  for (unsigned int i=1;i<splitchoices.size();i++) {
    if (splitchoices[i].reason_value>retval.reason_value)
      retval=splitchoices[i];
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

SplitChoice OneBoundIntervaluation::intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals) const
{
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
  Index tmp_idx;
  //cycle through the different parameter values for the given interval and determine branching objective values
  for (int i=0; i<npars;i++) {
    int_splitchoices[2*i].parameterID = parameterIDs[i];
    int_splitchoices[2*i].intervalID = intervalID;
    int_splitchoices[2*i+1].parameterID = parameterIDs[i];
    int_splitchoices[2*i+1].intervalID = intervalID;
    for (int j=0;j<intervals.size();j++) {
      if (ID_vec[j]==intervalID && para_vec[j]==parameterIDs[i]) {
	// parameter entry for interval in question found!
	tmp_idx = intervals.getOtherBndIdx(j);
	if (!dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(j)))->IsHomogeneous()) {
	  const Number* sensitivity_value_1 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(j)))->Values();
	  int_splitchoices[2*i].reason_value = *(sensitivity_value_1+controlrow)*intervalwidths[i];
	} else {
	  const Number sensitivity_value_1 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(j)))->Scalar();
	  int_splitchoices[2*i].reason_value = sensitivity_value_1*intervalwidths[i];
	}
	if (!dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(tmp_idx)))->IsHomogeneous()) {
	  const Number* sensitivity_value_2 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(tmp_idx)))->Values();
	  int_splitchoices[2*i+1].reason_value = *(sensitivity_value_2+controlrow)*intervalwidths[i];
	} else {
	  const Number sensitivity_value_2 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(tmp_idx)))->Scalar();
	  int_splitchoices[2*i+1].reason_value = sensitivity_value_2*intervalwidths[i];
	}
	j=intervals.size();
      }
    }
  }
  //let intbrancher decide, which parameter is critical on this interval
  retval = int_brancher->branchInterval(int_splitchoices);
  return retval;
}

SplitChoice BothBoundIntervaluation::intervaluateInterval(const Index& controlrow, const Index& intervalID,SmartPtr<MultiVectorMatrix> mv_sens,SmartPtr<OptionsList> options, const IntervalInfoSet& intervals) const
{
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
    for (int j=0;j<intervals.size();j++) {
      if (ID_vec[j]==intervalID && para_vec[j]==parameterIDs[i]) {
	// parameter entry for interval in question found!
	tmp_idx = intervals.getOtherBndIdx(j);
	if (!dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(j)))->IsHomogeneous()) {
	  const Number* sensitivity_value_1 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(j)))->Values();
	  tmp_value = *(sensitivity_value_1+controlrow);
	} else {
	  const Number sensitivity_value_1 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(j)))->Scalar();
	  tmp_value = sensitivity_value_1;
	}
	if (!dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(tmp_idx)))->IsHomogeneous()) {
	  const Number* sensitivity_value_2 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(tmp_idx)))->Values();
	  int_splitchoices[i].reason_value = tmp_value*(*(sensitivity_value_2+controlrow))*intervalwidths[i];
	} else {
	  const Number sensitivity_value_2 = dynamic_cast<const DenseVector*>(GetRawPtr(mv_sens->GetVector(tmp_idx)))->Scalar();
	  int_splitchoices[i].reason_value = tmp_value*(sensitivity_value_2)*intervalwidths[i];
	}
	j=intervals.size();
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
  SmartPtr<MultiVectorMatrix> mv_sens = dynamic_cast<MultiVectorMatrix*>(GetRawPtr(getSensitivityMatrix(app)));
  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const DenseVector> p = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  IntervalInfoSet intervals = IntervalInfoSet(p);

  std::vector<SplitChoice> splitchoices = branchmode->branchSensitivityMatrix(mv_sens,options,intervals);

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
    } if (sensemode == "GMRES" || sensemode == "MINRES") {
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
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));

  KKT_ = pd_search->PDSolver();
  jnl_ = app->Jnlst();

  //assign local original Ipopt data
  x_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->x()));
  s_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->s()));
  y_c_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->y_c()));
  y_d_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->y_d()));
  z_L_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->z_L()));
  z_U_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->z_U()));
  v_L_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->v_L()));
  v_U_ = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->v_U()));

  SmartPtr<const DenseVector>  p_ = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));

  Wp_ = orig_nlp->h_p(*x_, 1.0, *y_c_, *y_d_);
  Ap_  = orig_nlp->jac_c_p(*x_);
  Bp_  = orig_nlp->jac_d_p(*x_);
  W_ = orig_nlp->h(*x_, 1.0, *y_c_, *y_d_);
  A_ = orig_nlp->jac_c(*x_);
  B_ = orig_nlp->jac_d(*x_);

  x_L_  = orig_nlp->x_L();
  Px_L_ = orig_nlp->Px_L();
  x_U_  = orig_nlp->x_U();
  Px_U_ = orig_nlp->Px_U();
  d_L_  = orig_nlp->d_L();
  Pd_L_ = orig_nlp->Pd_L();
  d_U_  = orig_nlp->d_U();
  Pd_U_ = orig_nlp->Pd_U();

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
  n_u_ = int(u_indices_.size());

  // temporary for top, bot and rhs (=top+bot) dim determination
  const Index  n_x = int(x_intervalIDs_.size())-n_u_;
  const Index  n_s = int(d_intervalIDs_.size());
  const Index  n_c = int(c_intervalIDs_.size());
  const Index  n_d = int(d_intervalIDs_.size());
  const Index n_zl = Px_L_->NCols();
  const Index n_zu = Px_U_->NCols();
  const Index n_vl = Pd_L_->NCols();
  const Index n_vu = Pd_U_->NCols();

  top_dim_ = n_x + n_u_ + n_s + n_c + n_d + n_zl + n_zu + n_vl + n_vu;

  // determine number of upper and lower bounds on controls - must be removed from bot
  n_ul_ = 0;
  n_uu_ = 0;

  const Index* x_l_pos = dynamic_cast<const ExpansionMatrix*>(GetRawPtr(Px_L_))->ExpandedPosIndices();
  const Index* x_u_pos = dynamic_cast<const ExpansionMatrix*>(GetRawPtr(Px_U_))->ExpandedPosIndices();

  for (unsigned int i=0;i<u_indices_.size(); i++) {
    for (int j=0;j<Px_L_->NCols(); j++) {
      if (*(x_l_pos+j) == u_indices_[i])
	n_ul_++;
    }
    for (int j=0;j<Px_U_->NCols(); j++) {
      if (*(x_u_pos+j) == u_indices_[i])
	n_uu_++;
    }
  }

  const Index bot_dim = int((n_x + n_s + n_c + n_d + n_zl - n_ul_ + n_zu - n_uu_ + n_vu + n_vl) / n_i_);
  rhs_dim_ = top_dim_+bot_dim;

  // printf("\nPx_L_: NRows: %d  NCols: %d",Px_L_->NRows(),Px_L_->NCols());
  // printf("\nPx_U_: NRows: %d  NCols: %d",Px_U_->NRows(),Px_U_->NCols());
  // printf("\nPd_L_: NRows: %d  NCols: %d",Pd_L_->NRows(),Pd_L_->NCols());
  // printf("\nPd_U_: NRows: %d  NCols: %d",Pd_U_->NRows(),Pd_U_->NCols());
  // printf("\nz_L_: Dim: %d",z_L_->Dim());
  // printf("\nz_U_: Dim: %d",z_U_->Dim());
  // printf("\nv_L_: Dim: %d",v_L_->Dim());
  // printf("\nv_U_: Dim: %d",v_U_->Dim());
  // printf("\nx_L_: Dim: %d",x_L_->Dim());
  // printf("\nx_U_: Dim: %d",x_U_->Dim());
  // printf("\nd_L_: Dim: %d",d_L_->Dim());
  // printf("\nd_U_: Dim: %d",d_U_->Dim());

  // parts for the symmetricized MINRES

  // set up slack values for x bounds
  Sl_x_L_ = x_->MakeNewDenseVector();
  Px_L_->MultVector(1.0,*x_L_,0.0,*Sl_x_L_);

  Sl_x_U_ = x_->MakeNewDenseVector();
  Px_U_->MultVector(1.0,*x_U_,0.0,*Sl_x_U_);

  Sl_x_L_->AddOneVector(1.0,*x_,-1.0);
  Sl_x_U_->AddOneVector(-1.0,*x_,1.0);

  // printf("\n");
  // Sl_x_L_->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"Slxl");
  // printf("\n");
  // Sl_x_U_->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"Slxu");

  // expand bound multipliers for x
  z_L_DIV_Sl_ = x_->MakeNewDenseVector();
  Px_L_->MultVector(1.0,*z_L_,0.0,*z_L_DIV_Sl_);

  z_U_DIV_Sl_ = x_->MakeNewDenseVector();
  Px_U_->MultVector(1.0,*z_U_,0.0,*z_U_DIV_Sl_);

  z_L_DIV_Sl_->ElementWiseDivide(*Sl_x_L_);
  z_U_DIV_Sl_->ElementWiseDivide(*Sl_x_U_);

  // set up slack values for inequality bounds
  Sl_s_L_ = s_->MakeNewDenseVector();
  Pd_L_->MultVector(1.0,*d_L_,0.0,*Sl_s_L_);

  Sl_s_U_ = s_->MakeNewDenseVector();
  Pd_U_->MultVector(1.0,*d_U_,0.0,*Sl_s_U_);

  Sl_s_L_->AddOneVector(1.0,*s_,-1.0);
  Sl_s_U_->AddOneVector(-1.0,*s_,1.0);

  printf("\n");
  Sl_s_L_->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"Slsl");
  printf("\n");
  Sl_s_U_->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"Slsu");

  // expand bound multipliers for inequalities
  v_L_DIV_Sl_ = s_->MakeNewDenseVector();
  Pd_L_->MultVector(1.0,*v_L_,0.0,*v_L_DIV_Sl_);

  v_U_DIV_Sl_ = s_->MakeNewDenseVector();
  Pd_U_->MultVector(1.0,*v_U_,0.0,*v_U_DIV_Sl_);

  v_L_DIV_Sl_->ElementWiseDivide(*Sl_s_L_);
  v_U_DIV_Sl_->ElementWiseDivide(*Sl_s_U_);

  printf("\n");
  z_L_DIV_Sl_->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"zldiv");
  printf("\n");
  z_U_DIV_Sl_->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"zudiv");
  printf("\n");
  v_L_DIV_Sl_->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"vldiv");
  printf("\n");
  v_U_DIV_Sl_->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"vudiv");
}

/* apply GMRES-approximation of a split to all intervals and decide for the smartest split*/
SplitDecision LinearizeKKT::applySplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  // testGMRES();
  SmartPtr<OptionsList> options = app->Options();
  SplitDecision retval;

  // provde MetaData for sens_vspace (needed by branchSensitivityMatrix())
  std::vector<Index> int_IDs = collectIntMetaData();

  SmartPtr<DenseVectorSpace> sens_vspace = new DenseVectorSpace(rhs_dim_);
  sens_vspace->SetIntegerMetaData("intervalID",int_IDs);
  SmartPtr<DenseVector> sens_vec = sens_vspace->MakeNewDenseVector();
  SmartPtr<MultiVectorMatrixSpace> sens_space = new MultiVectorMatrixSpace(intervals_.size(),*sens_vspace);
  SmartPtr<MultiVectorMatrix> mv_sens = sens_space->MakeNewMultiVectorMatrix();

  std::string sensemode;
  if (options->GetStringValue("sensemode",sensemode ,"")){

    for (int i=0;i<intervals_.size();i++) {
      // get GMRES approximated split results wrt each single interval
      // i+1: in case intervalIDs start with 1 (which they do)
      if (intervals_.isUpper(i)) {

	if (sensemode =="MINRES")
	  sens_vec = applyMINRESOnInterval(app,intervals_.getIntervalID(i),intervals_.getParameterID(i));
	else
	  sens_vec = applyGMRESOnInterval(app,intervals_.getIntervalID(i),intervals_.getParameterID(i));
	//sens_vec->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"sense_vec");
	// printf("\n\n");
      } else
	sens_vec->Set(0);

      mv_sens->SetVector(i,*sens_vec);
      sens_vec = sens_vspace->MakeNewDenseVector();
    }
  // printf("\n\n");
  // mv_sens->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"mv_sens");

  // cycle through var space interval flags to identify and save control indexes
  std::vector<Index> ctrl_rows;
  assert(dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->HasIntegerMetaData("intervalID"));
  const std::vector<Index> var_int_flags = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(mv_sens->ColVectorSpace()))->GetIntegerMetaData("intervalID");
  for (int i=0;i<x_->Dim();i++) {
    if (!var_int_flags[i])
      ctrl_rows.push_back(i);
  }

  // chose the split with best results
  BranchingCriterion* branchmode = assignBranchingMethod(options);
  std::vector<SplitChoice> splitchoices = branchmode->branchSensitivityMatrix(mv_sens,options,intervals_,ctrl_rows,true);
  ControlSelector* pickfirst = assignControlMethod(options);
  retval = pickfirst->decideSplitControl(splitchoices);

  return retval;
  }

  assert(0);
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
  const unsigned int small_dim = int(indices.size());

  assert((unsigned) large_dim>=small_dim);

  if ((unsigned) large_dim==small_dim) {
//    printf("\nLinearizeKKT::expandVector(): WARNING: called to expand a vector to original size.");
    return original;
  } else {
    Index* exppos = new Index[small_dim];
    for (unsigned int i=0;i<small_dim;i++)
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
  const unsigned int large_dim = original->Dim();
  const unsigned int small_dim = int(indices.size());

  assert(large_dim>=small_dim);
  if (large_dim==small_dim) {
    //    printf("\nLinearizeKKT::shrink(): WARNING: called to shrink a vector to original size.");
    return original;
  } else {

    Index* exppos = new Index[small_dim];
    for (unsigned int i=0;i<small_dim;i++)
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

SmartPtr<MultiVectorMatrix> LinearizeKKT::computeVMultEye(SmartPtr<const DenseVector> target) const
{
  // prepare spaces and variables
  SmartPtr<DenseVectorSpace> retvec_space = new DenseVectorSpace(target->Dim());
  SmartPtr<DenseVector> retvec;
  SmartPtr<MultiVectorMatrixSpace> retval_space = new MultiVectorMatrixSpace(target->Dim(),*retvec_space);
  SmartPtr<MultiVectorMatrix> retval = retval_space->MakeNewMultiVectorMatrix();
  const bool homog = target->IsHomogeneous();
  Number val;
  Number* vals = new Number[target->Dim()];
  Number* retvals = new Number[target->Dim()];

  // get target values
  if (homog) {
    val = target->Scalar();
  } else {
    const Number* tmp_vals = target->Values();
    std::copy(tmp_vals,tmp_vals+target->Dim(),vals);
  }

  // set up vector values and assign vectors
  for (int i=0;i<target->Dim();i++) {
    retvec = retvec_space->MakeNewDenseVector();
    for (int j=0;j<target->Dim();j++)
      retvals[j] = 0;

    if (homog)
      retvals[i] = val;
    else
      retvals[i] = vals[i];

    retvec->SetValues(retvals);
    retval->SetVector(i,*retvec);
  }

  return retval;
}

SmartPtr<const DenseVector> LinearizeKKT::computeDeltaP(const IntervalInfoSet& intervals,const Index& interval, const Index& parameterID) const
{
  // collect neccessary data
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
  Number int_width;

  // determine interval width for given parameter
  for (unsigned int j=0; j<parameterIDs.size(); j++) {
    if (parameterIDs[j] == parameterID) {
      int_width = int_widths[j];
      break;
    }
  }

  // for (int i=0;i<int_widths.size();i++)
  //   printf("\nint_widths[%d] = %f",i,int_widths[i]);

  // push back desired values for the given parameterID, other than that 0
  for (unsigned int k=0; k<intervalflags.size(); k++) {
    if (intervalflags[k] == interval) {
      if (parameterflags[k] == parameterID) {
	if (intervals.isUpper(k))
	  ret_values[counter] = -0.5*int_width;
	else
	  ret_values[counter] = 0.5*int_width;
	counter++;
      } else {
      ret_values[counter] = 0;
      counter++;
      }
    }
  }

  retval->SetValues(ret_values);
  return dynamic_cast<const DenseVector*>(GetRawPtr(retval));
}

SmartPtr<DenseVector> LinearizeKKT::applyGMRESOnInterval(SmartPtr<IpoptApplication> app,const Index& interval, const Index& parameter)
{  /*
  // initialize neccessary quantities
  const Index n_parameters = intervals_.getParameterCount();
  std::vector<Index> intervalIDs = intervals_.getIntervalIDVec();
  SmartPtr<DenseVectorSpace> z_space = new DenseVectorSpace(rhs_dim_);
  SmartPtr<DenseVector> z_postsplit = z_space->MakeNewDenseVector();
  SmartPtr<const DenseVector> deltap;

  SmartPtr<DenseVectorSpace> res_vspace = new DenseVectorSpace(rhs_dim_);
  SmartPtr<DenseVector> res_vsense = res_vspace->MakeNewDenseVector();

  SmartPtr<const IteratesVectorSpace> top_space;
  SmartPtr<IteratesVector> rhs_top;
  SmartPtr<DenseVector> rhs_x,rhs_y_c,rhs_y_d;
  SmartPtr<ShiftVector> rhs;
  SmartPtr<ShiftVector> z_cond;
  SmartPtr<ShiftVector> lhs;
  SmartPtr<ShiftVector> x0;
  Number tolo;
  Index n_m,n_rst;

  // init sensitivity matrix and -space
  SmartPtr<MultiVectorMatrixSpace> sense_space = new MultiVectorMatrixSpace(n_parameters*2,*z_space);
  SmartPtr<MultiVectorMatrix> sense = sense_space->MakeNewMultiVectorMatrix();
  SmartPtr<MultiVectorMatrix> res_sense = sense_space->MakeNewMultiVectorMatrix();

  // start algorithm on each existing intervalID
  Index s_count = 0;
  for (unsigned  int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i]==interval) {

      // split interval indices into to be shifted interval entries and remainder
      if (!splitIntervalIndices(x_intervalIDs_,shift_x_indices_,interval))
	printf("\nLinearizeKKT::applyGMRESOnInterval(): ERROR: unable to split x_intervalIDs_");
      if (!splitIntervalIndices(c_intervalIDs_,shift_c_indices_,interval))
	printf("\nLinearizeKKT::applyGMRESOnInterval(): ERROR: unable to split c_intervalIDs_");
      if (!splitIntervalIndices(d_intervalIDs_,shift_d_indices_,interval))
	printf("\nLinearizeKKT::applyGMRESOnInterval(): ERROR: unable to split d_intervalIDs_");

      // get different rhs entries for this column
      x_i_ = extractColumn(Wp_,i);
      y_c_i_ = extractColumn(Ap_,i);
      y_d_i_ = extractColumn(Bp_,i);

      // assign rhs subvectors
      u_i_ = shrink(x_i_,u_indices_);

      // create ShiftVector for rhs
      top_space = dynamic_cast<const IteratesVectorSpace*>(GetRawPtr(app->IpoptDataObject()->curr()->OwnerSpace()));
      rhs_top = top_space->MakeNewIteratesVector();
      rhs_top->Set(0.0);
      rhs_top->Set_x_NonConst(*x_i_);
      rhs_top->Set_y_c_NonConst(*y_c_i_);
      rhs_top->Set_y_d_NonConst(*y_d_i_);

      SmartPtr<DenseVectorSpace> x_sp = new DenseVectorSpace(shift_x_indices_.size());
      rhs_x = x_sp->MakeNewDenseVector();
      rhs_x = shrink(x_i_,shift_x_indices_);
      SmartPtr<DenseVectorSpace> y_c_sp = new DenseVectorSpace(shift_c_indices_.size());
      rhs_y_c = y_c_sp->MakeNewDenseVector();
      rhs_y_c = shrink(y_c_i_,shift_c_indices_);
      SmartPtr<DenseVectorSpace> y_d_sp = new DenseVectorSpace(shift_d_indices_.size());
      rhs_y_d = y_d_sp->MakeNewDenseVector();
      rhs_y_d = shrink(y_d_i_,shift_d_indices_);

      rhs = new ShiftVector(rhs_top,rhs_x,rhs_y_c,rhs_y_d);

      // create x0
      x0 = new ShiftVector(*rhs);
      x0->Set(0.0);

      // compute solution of the GMRES-Step
      assignGMRESOptions(app->Options(),tolo,n_m,n_rst);
      z_cond = computeGMRES(rhs,x0,tolo,n_m);

      // undo conditionning step to get the approximation of the postsplit sensitivity column
      lhs = computePMultVector(z_cond);
      z_postsplit = lhs->getDVector();

      // add the resulting sense vector to matrix ....
      sense->SetVector(s_count,*z_postsplit);
      s_count++;
    }
  }

  // get the parameter perturbation for this split
  deltap = computeDeltaP(intervals_,interval,parameter);
  // calculate absolute sensitivity values
  sense->MultVector(1.0,*deltap,0.0,*res_vsense);

  return res_vsense; */
}

/** given an Index-type vector containing a list of intervalIDs, loads indices matching the given interval into the shift_indices part, the rest into static_indices **/
bool LinearizeKKT::splitIntervalIndices(const std::vector<Index>& intervalIDs,std::vector<Index>& shift_indices,const Index& interval)
{
  bool retval = false;
  // make sure nothings in the way of push_backs
  shift_indices.clear();

  assert(intervalIDs.size());

  // cycle through intervalIDs and assign the indices
  for (unsigned int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i] == interval) {
      shift_indices.push_back(i);
    }
  }

  if (shift_indices.size()>0)
    retval = true;
  else {
    printf("\nLinearizeKKT::splitIntervalIndices(): looking for intervalID %d in the following list:",interval);
  for (unsigned int i=0;i<intervalIDs.size();i++)
    printf("\n  intervalIDs(%d) = %d",i,intervalIDs[i]);
  printf("\nwas not successful.");
  }

  return retval;
}

void LinearizeKKT::assignGMRESOptions(SmartPtr<OptionsList> options,Number& tol, Index& n_max, Index& n_rest) const
{
  if (options->GetNumericValue("gmr_tol",tol,"")) {
    // printf("\nLinearizeKKT::assignGMRESOptions(): gmr_tol set to: %e",tol);
  } else {
    //    printf("\nLinearizeKKT::assignGMRESOptions(): gmr_tol set to default value");
    tol = 1.0e-5;
  }
  if (options->GetIntegerValue("gmr_n_max",n_max,"")) {
    //    printf("\nLinearizeKKT::assignGMRESOptions(): n_max set to: %d",n_max);
  } else {
    //    printf("\nLinearizeKKT::assignGMRESOptions(): n_max set to default value");
    n_max = 20;
  }

  // restart not implemented yet

 /*  if (options->GetIntegerValue("gmr_n_rst",n_rest,"")) {
    printf("\nLinearizeKKT::assignGMRESOptions(): n_rest set to: %d",n_rest);
  } else {
    printf("\nLinearizeKKT::assignGMRESOptions(): n_rest set to default value");
    n_rest = 0;
    } */
}

SmartPtr<ShiftVector> LinearizeKKT::computeGMRES(SmartPtr<ShiftVector>b,SmartPtr<ShiftVector>x0,const Number& tol, const Index& n_max, const Index& n_rest)
{
  // initialize neccessary vars
  SmartPtr<ShiftVector> r = computeR(b,x0);
  SmartPtr<ShiftVector> x = new ShiftVector(*x0);

  std::vector<SmartPtr<ShiftVector> > v(n_max+1);
  std::vector<SmartPtr<ShiftVector> > w(n_max+1);

  Number beta;
  Index g_cnt = 0;
  std::vector<Number> h((n_max+1)*(n_max+1));
  std::vector<Number> s(n_max+1);
  std::vector<Number> c(n_max+1);
  std::vector<Number> gamma(n_max+1);
  std::vector<Number> y(n_max);

  SmartPtr<ShiftVector> tmp;

  // for residual output
  std::string fname = "residuals.dat";
  std::ofstream residuals;
  residuals.open(fname.c_str());
  char buffer[63];
  residuals << "#.dat file automatically generated by AMPL intervallization routine\n#Ben Waldecker Sep 2012\n";

  // preloop inital var values
  v[0]= computeV(r);
  gamma[0] = r->Nrm2();

  // outer loop
  for (Index j=0; j<n_max;j++) {

    w[j] = computeKPIMultVector(v[j]);

    // update h[i,j] entries
    for (Index i=0; i<=j;i++) {
      h[i+j*n_max] = v[i]->Dot(*w[j]);
    }

    // calculate w[j]
    for (Index i=0;i<=j;i++) {
      w[j]->AddOneVector(-1.0*h[i+j*n_max],*v[i],1.0);
    }

    // rotate it
    for (Index i=0;i<j;i++) {
      double hij = h[i+j*n_max];
      double hij1 = h[i+1+j*n_max];
      h[i+j*n_max] = c[i+1]*hij + s[i+1]*hij1;
      h[i+1+j*n_max] = -s[i+1]*hij + c[i+1]*hij1;
    }

    // expand h, first step
    h[j+1+j*n_max] = w[j]->Nrm2();

    // expand h, second step -
    beta = sqrt(h[j+j*n_max]*h[j+j*n_max]+h[j+1+j*n_max]*h[j+1+j*n_max]);
    s[j+1] = h[j+1+j*n_max]/beta;
    c[j+1] = h[j+j*n_max]/beta;
    h[j+j*n_max] = beta;

    // expand and update gammas
    gamma[j+1] = -s[j+1]*gamma[j];
    gamma[j] = c[j+1]*gamma[j];
    g_cnt++;

    // check quality of current solution
    if ( (fabs(gamma[j+1]) > tol) && ( !n_max || (n_max && (j!=n_max-1)) ) ) {
      // insufficient - update v[j+1]
      assert(h[j+1+j*n_max]);
      v[j+1] = new ShiftVector(*w[j]);
      v[j+1]->Scal(1/h[j+1+j*n_max]);

      // target accuracy or maximum number of iterations achieved - calculate solution
    } else {
      // set up multipliers y[i]
      for (Index i=j;i>=0;i--) { // check this!
	assert(h[i+i*n_max]);

	// set up sum[k] h[i+k*n_max] y[k]
	beta = 0;
	for (Index k=i+1;k<=j;k++) {
	  beta += h[i+k*n_max]*y[k];
	}
	y[i]= (1/h[i+i*n_max])*(gamma[i]-beta);
      }

      // compute solution
      for (Index i=0; i<=j;i++) {
	x->AddTwoVectors(y[i],*v[i],0.0,*v[i],1.0);
      }

      break;
    }

    /* set up current solution to calculate absolute residuum and write it to file
       in order to check for convergence */
    for (Index i=j;i>=0;i--) {
      assert(h[i+i*n_max]);
      beta = 0;
      for (Index k=i+1;k<=j;k++) {
	beta += h[i+k*n_max]*y[k];
      }
      y[i]= 1/h[i+i*n_max]*(gamma[i]-beta);
    }
    for (Index i=0; i<=j;i++) {
      x->AddTwoVectors(y[i],*v[i],0.0,*v[i],1.0);
    }
    SmartPtr<ShiftVector> fooo = computeKPIMultVector(x);
    fooo->AddOneVector(-1.0,*b,1.0);
    Number residual = fooo->Nrm2();

    //write residual into residuals.dat file for python access
    sprintf(buffer,"\nresidual[%d] = %e,      (relative: %e)",j,residual, residual/b->Nrm2());
    residuals << buffer;
    x = new ShiftVector(*x0);
  }
  residuals << "\n\n#end of file";
  residuals.close();

  //write gamma[j] into gamma.dat file for python access
  fname = "gammas.dat";
  std::ofstream gammas;
  gammas.open(fname.c_str());
  gammas << "#.dat file automatically generated by AMPL intervallization routine\n#Ben Waldecker Sep 2012\n";
  for (int i=0;i<g_cnt;i++) {
    sprintf(buffer,"\ngamma[%d] = %e",i,gamma[i]);
    gammas << buffer;
  }
  gammas << "\n\n#end of file";
  gammas.close();

  return x;
}

SmartPtr<ShiftVector> LinearizeKKT::computeR(SmartPtr<ShiftVector> b, SmartPtr<ShiftVector>x) const
{
  SmartPtr<ShiftVector> retval = new ShiftVector(*b);
  SmartPtr<ShiftVector> Ax = computeKPIMultVector(x);
  retval->AddOneVector(-1.0,*Ax,1.0);

  return retval;
}

SmartPtr<ShiftVector> LinearizeKKT::computeV(SmartPtr<ShiftVector> r) const
{
  SmartPtr<ShiftVector> retval = new ShiftVector(*r);

  const Number rr = r->Nrm2();
  assert(rr);
  retval->Scal(1/rr);

  return retval;
}

SmartPtr<ShiftVector> LinearizeKKT::computeAMultVector(SmartPtr<ShiftVector> target) const
{
  // hijacked by MINRES
}

SmartPtr<ShiftVector> LinearizeKKT::computeBMultVector(SmartPtr<ShiftVector> target) const
{
  // hijacked by MINRES
}

SmartPtr<ShiftVector> LinearizeKKT::computeCMultVector(SmartPtr<ShiftVector> target) const
{
  // hijacked by MINRES
}

SmartPtr<ShiftVector> LinearizeKKT::computeDMultVector(SmartPtr<ShiftVector> target) const
{
  // hijacked by MINRES
}

SmartPtr<ShiftVector> LinearizeKKT::computeAMVminres(SmartPtr<ShiftVector> target) const
{
  // set up containers for Matrix Vector Multiplication and retval part,
  SmartPtr<DenseVector> retval_x = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  SmartPtr<DenseVector> extractor = retval_x->MakeNewDenseVector();
  SmartPtr<DenseVector> retval_s = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  SmartPtr<DenseVector> retval_y_c = dynamic_cast<DenseVector*>(target->top()->y_c()->MakeNew());
  SmartPtr<DenseVector> retval_y_d = dynamic_cast<DenseVector*>(target->top()->y_d()->MakeNew());

  // first large KKT-matrix row
  W_->MultVector(1.0,*target->top()->x(),0.0,*extractor);
  retval_x->AddOneVector(1.0,*extractor,0.0);

  // // additional part due to symmetricising of the KKT matrix
  // extractor = retval_x->MakeNewDenseVector();
  // extractor->AddTwoVectors(1.0,*z_L_DIV_Sl_,1.0,*z_U_DIV_Sl_,0.0);
  // extractor->ElementWiseMultiply(*target->top()->x());
  // retval_x->AddOneVector(1.0,*extractor,1.0);

  // extractor = retval_x->MakeNewDenseVector();
  // A_->TransMultVector(1.0,*target->top()->y_c(),0.0,*extractor);
  // retval_x->AddOneVector(1.0,*extractor,1.0);

  // extractor = retval_x->MakeNewDenseVector();
  // B_->TransMultVector(1.0,*target->top()->y_d(),0.0,*extractor);
  // retval_x->AddOneVector(1.0,*extractor,1.0);

  // // second large KKT-matrix row
  // extractor = retval_s->MakeNewDenseVector();
  // extractor->AddTwoVectors(1.0,*v_L_DIV_Sl_,1.0,*v_U_DIV_Sl_,0.0);
  // extractor->ElementWiseMultiply(*target->top()->s());
  // retval_s->AddTwoVectors(-1.0,*target->top()->y_d(),1.0,*extractor,0.0);

  // third large KKT-matrix row
  extractor = retval_y_c->MakeNewDenseVector();
  A_->MultVector(1.0,*target->top()->x(),0.0,*extractor);
  retval_y_c->AddOneVector(1.0,*extractor,0.0);

  // fourth large KKT-matrix row
  extractor = retval_y_d->MakeNewDenseVector();
  B_->MultVector(1.0,*target->top()->x(),0.0,*extractor);
  retval_y_d->AddTwoVectors(1.0,*extractor,-1.0,*target->top()->s(),0.0);

  //set up retval
  SmartPtr<ShiftVector> retval = new ShiftVector(*target);
  retval->Scal(0.0);
  retval->Set_x_top(*retval_x);
  //  retval->Set_s_top(*retval_s);
  retval->Set_y_c_top(*retval_y_c);
  retval->Set_y_d_top(*retval_y_d);

  return new ShiftVector(*retval);
}

SmartPtr<ShiftVector> LinearizeKKT::computeBMVminres(SmartPtr<ShiftVector> target) const
{
  SmartPtr<DenseVectorSpace> u_space = new DenseVectorSpace(int(u_indices_.size()));
  SmartPtr<DenseVector> retval_u = u_space->MakeNewDenseVector();

  // get (W_i,shift times rhs_i) with i= (shift_h, u)
  SmartPtr<DenseVector> tmp_rhs = expand(target->x(),shift_x_indices_,W_->NCols());
  SmartPtr<DenseVectorSpace> extractor_space = new DenseVectorSpace(W_->NRows());
  SmartPtr<DenseVector> extractor = extractor_space->MakeNewDenseVector();
  W_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map results back to get u part
  SmartPtr<DenseVector> extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,0.0);

  // additional part due to symmetricising of the KKT matrix
  extractor = extractor_space->MakeNewDenseVector();
  extractor->AddTwoVectors(1.0,*z_L_DIV_Sl_,1.0,*z_U_DIV_Sl_,0.0);
  extractor->ElementWiseMultiply(*target->top()->x());
  //map results back to get u part
  extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,1.0);

  //get (A_shift^T times rhs_i) parts with i= (u, shift_c)
  tmp_rhs = expand(target->y_c(),shift_c_indices_,A_->NRows());
  extractor_space = new DenseVectorSpace(A_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  A_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map results back to get u part
  extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,1.0);

  //get (B_shift^T times rhs_i) parts with i= (u, shift_d)
  tmp_rhs = expand(target->y_d(),shift_d_indices_,B_->NRows());
  extractor_space = new DenseVectorSpace(B_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  B_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map results back to get u part
  extr_u_part = shrink(extractor,u_indices_);
  retval_u->AddOneVector(1.0,*extr_u_part,1.0);

  SmartPtr<ShiftVector> retval = new ShiftVector(*target);
  retval->Scal(0.0);
  // insert calculated u part into retval
  SmartPtr<DenseVector> top_x = expand(retval_u,u_indices_,target->top()->x()->Dim());
  retval->Set_x_top(*top_x);

  return new ShiftVector(*retval);
}

SmartPtr<ShiftVector> LinearizeKKT::computeCMVminres(SmartPtr<ShiftVector> target) const
{
  SmartPtr<DenseVector> retval_x = dynamic_cast<DenseVector*>(target->top()->x()->MakeNewCopy());
  SmartPtr<DenseVector> retval_u = shrink(retval_x,u_indices_);
  retval_x = expand(retval_u,u_indices_,W_->NCols());

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
  SmartPtr<DenseVector> retval_y_c = shrink(retval_c_nonconst,shift_c_indices_);
  SmartPtr<DenseVector> retval_y_d = shrink(retval_d_nonconst,shift_d_indices_);

  //setup return value
  SmartPtr<ShiftVector> retval = new ShiftVector(*target);
  retval->Scal(0.0);
  retval->Set_x(*retval_x);
  retval->Set_y_c(*retval_y_c);
  retval->Set_y_d(*retval_y_d);

  return new ShiftVector(*retval);
}

SmartPtr<ShiftVector> LinearizeKKT::computeDMVminres(SmartPtr<ShiftVector> target) const
{
  // init retval parts
  SmartPtr<DenseVector> retval_x = target->x()->MakeNewDenseVector();
  SmartPtr<DenseVector> retval_s = dynamic_cast<DenseVector*>(target->s()->MakeNewCopy());
  SmartPtr<DenseVector> retval_y_c = target->y_c()->MakeNewDenseVector();
  SmartPtr<DenseVector> retval_y_d = target->y_d()->MakeNewDenseVector();

  // first row
  // get W_shift,shift times rhs_shift_h part
  SmartPtr<DenseVector> tmp_rhs = expand(target->x(),shift_x_indices_,W_->NCols());
  SmartPtr<DenseVectorSpace> extractor_space = new DenseVectorSpace(W_->NRows());
  SmartPtr<DenseVector> extractor = extractor_space->MakeNewDenseVector();
  W_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  // map results back
  SmartPtr<DenseVector> extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,0.0);

  // additional part due to symmetricising of the KKT matrix
  extractor = extractor_space->MakeNewDenseVector();
  extractor->AddTwoVectors(1.0,*z_L_DIV_Sl_,1.0,*z_U_DIV_Sl_,0.0);
  extractor->ElementWiseMultiply(*target->top()->x());
  tmp_rhs = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extractor,1.0);

  // get A_shift,shift_c^T times rhs_shift_c part
  tmp_rhs = expand(target->y_c(),shift_c_indices_,A_->NRows());
  extractor_space = new DenseVectorSpace(A_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  A_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  // map results back
  extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,1.0);

  // get (B_shift^T times rhs_i) parts with i= (u, shift_d)
  tmp_rhs = expand(target->y_d(),shift_d_indices_,B_->NRows());
  extractor_space = new DenseVectorSpace(B_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  B_->TransMultVector(1.0,*tmp_rhs,0.0,*extractor);
  // map result back
  extr_shift_part = shrink(extractor,shift_x_indices_);
  retval_x->AddOneVector(1.0,*extr_shift_part,1.0);

  // second row
  extractor_space = new DenseVectorSpace(B_->NCols());
  extractor = extractor_space->MakeNewDenseVector();
  extractor->AddTwoVectors(1.0,*v_L_DIV_Sl_,1.0,*v_U_DIV_Sl_,0.0);
  extractor->ElementWiseMultiply(*target->top()->s());
  tmp_rhs = dynamic_cast<DenseVector*>(target->top()->y_d()->MakeNewCopy());
  retval_s = shrink(tmp_rhs,shift_d_indices_);
  tmp_rhs = shrink(extractor,shift_d_indices_);
  retval_s->AddOneVector(1.0,*tmp_rhs,-1.0);

  // third row
  // get (A_shift times rhs_shift_c) part
  tmp_rhs = expand(target->x(),shift_x_indices_,A_->NCols());
  extractor = tmp_rhs->MakeNewDenseVector();
  A_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  //map results back
  retval_y_c = shrink(extractor,shift_c_indices_);

  // fourth row
  // get (B_shift times rhs_shift_d) part
  tmp_rhs = expand(target->x(),shift_x_indices_,B_->NCols());
  extractor = tmp_rhs->MakeNewDenseVector();
  B_->MultVector(1.0,*tmp_rhs,0.0,*extractor);
  // map results back
  retval_y_d = shrink(extractor,shift_d_indices_);

  SmartPtr<ShiftVector> retval = new ShiftVector(*target);
  retval->Scal(0.0);
  retval->Set_x(*retval_x);
  retval->Set_s(*retval_s);
  retval->Set_y_c(*retval_y_c);
  retval->Set_y_d(*retval_y_d);

  return new ShiftVector(*retval);
}

SmartPtr<ShiftVector> LinearizeKKT::computeKSymMultVector(SmartPtr<ShiftVector> target) const
{
  SmartPtr<ShiftVector> retval = new ShiftVector(*target);

  SmartPtr<ShiftVector> AM = computeAMVminres(target);
  SmartPtr<ShiftVector> BM = computeBMVminres(target);
  SmartPtr<ShiftVector> CM = computeCMVminres(target);
  SmartPtr<ShiftVector> DM = computeDMVminres(target);

  // printf("\n");
  // AM->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"AM");
  // printf("\n");
  // BM->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"BM");
  // printf("\n");
  // CM->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"CM");
  // printf("\n");
  // DM->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"DM");

  retval->AddTwoVectors(1.0,*AM,1.0,*BM,0.0);
  retval->AddTwoVectors(1.0,*CM,1.0,*DM,1.0);

  return new ShiftVector(*retval);
}

SmartPtr<ShiftVector> LinearizeKKT::computeKIMultVector(SmartPtr<ShiftVector> target) const
{
  // init neccessary quantities
  SmartPtr<const IteratesVectorSpace> pr_space = dynamic_cast<const IteratesVectorSpace*>(GetRawPtr(target->top()->OwnerSpace()));
  SmartPtr<IteratesVector> preretval = pr_space->MakeNewIteratesVector();
  SmartPtr<ShiftVector> retval = new ShiftVector(*target);
  retval->Scal(0.0);

  // set up and solve
  KKT_->Solve(1.0, 0.0, *target->top(), *preretval);

  retval->Set_x_top(*dynamic_cast<const DenseVector*>(GetRawPtr(preretval->x())));
  retval->Set_y_c_top(*dynamic_cast<const DenseVector*>(GetRawPtr(preretval->y_c())));
  retval->Set_y_d_top(*dynamic_cast<const DenseVector*>(GetRawPtr(preretval->y_d())));

  return new ShiftVector(*retval);
}

SmartPtr<ShiftVector> LinearizeKKT::computeKPIMultVector(SmartPtr<ShiftVector> target) const
{
  SmartPtr<ShiftVector> KIM = computeKIMultVector(target);
  SmartPtr<ShiftVector> CKIM = computeCMultVector(KIM);

  SmartPtr<ShiftVector> BM = computeBMultVector(target);
  SmartPtr<ShiftVector> DM = computeDMultVector(target);

  SmartPtr<ShiftVector> IM = new ShiftVector(*target);
  IM->Scal(0.0);
  IM->Set_x_top(*dynamic_cast<const DenseVector*>(GetRawPtr(target->top()->x())));
  IM->Set_y_c_top(*dynamic_cast<const DenseVector*>(GetRawPtr(target->top()->y_c())));
  IM->Set_y_d_top(*dynamic_cast<const DenseVector*>(GetRawPtr(target->top()->y_d())));

  SmartPtr<ShiftVector> retval = new ShiftVector(*target);
  retval->AddTwoVectors(1.0,*BM,1.0,*CKIM,0.0);
  retval->AddTwoVectors(1.0,*IM,1.0,*DM,1.0);

  return new ShiftVector(*retval);
}

SmartPtr<ShiftVector> LinearizeKKT::computePMultVector(SmartPtr<ShiftVector> target) const
{
  SmartPtr<IteratesVector> retval_top = target->top()->MakeNewIteratesVectorCopy();

  // do actual backsolve
  KKT_->Solve(1.0, 0.0, *target->top(), *retval_top);
  SmartPtr<const Vector> retval_x = retval_top->x()->MakeNewCopy();
  SmartPtr<const Vector> retval_c = retval_top->y_c()->MakeNewCopy();
  SmartPtr<const Vector> retval_d = retval_top->y_d()->MakeNewCopy();

  retval_top->Set(0.0);
  retval_top->Set_x(*retval_x);
  retval_top->Set_y_c(*retval_c);
  retval_top->Set_y_d(*retval_d);

  // set up retval

  SmartPtr<DenseVector> x = dynamic_cast<DenseVector*>(target->x()->MakeNewCopy());
  SmartPtr<DenseVector> s = dynamic_cast<DenseVector*>(target->s()->MakeNewCopy());
  SmartPtr<DenseVector> y_c = dynamic_cast<DenseVector*>(target->y_c()->MakeNewCopy());
  SmartPtr<DenseVector> y_d = dynamic_cast<DenseVector*>(target->y_d()->MakeNewCopy());
  SmartPtr<DenseVector> z_L = dynamic_cast<DenseVector*>(target->z_L()->MakeNewCopy());
  SmartPtr<DenseVector> z_U = dynamic_cast<DenseVector*>(target->z_U()->MakeNewCopy());
  SmartPtr<DenseVector> v_L = dynamic_cast<DenseVector*>(target->v_L()->MakeNewCopy());
  SmartPtr<DenseVector> v_U = dynamic_cast<DenseVector*>(target->v_U()->MakeNewCopy());

  SmartPtr<ShiftVector> retval = new ShiftVector(retval_top,x,s,y_c,y_d,z_L,z_U,v_L,v_U);

  return retval;
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

    /*
    SplitAlgorithm* splatter = new SplitIntAtBound();
    SplitDecision resulting_split = splatter->applySplitAlgorithm(app);
   */




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



SplitDecision SplitIntAtBound::applySplitAlgorithm(SmartPtr<IpoptApplication> app)
{
  SmartPtr<const DenseVector> x = dynamic_cast<const DenseVector*>(GetRawPtr(app->IpoptDataObject()->curr()->x()));

  SmartPtr<const DenseVectorSpace> x_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x->OwnerSpace()));

  assert(x_space->HasStringMetaData("idx_names"));
  const std::vector<std::basic_string<char> > names = x_space->GetStringMetaData("idx_names");

  std::vector<std::string> bound_values;
  std::vector<Index> indices;

  for (Index i=0;i<x->Dim();i++) {
    printf("\nname[%d]: %s comparision value: %d",i,names[i].c_str(),names[i].compare("xCU["));
    if (names[i].compare("xCU[")>0) {
      bound_values.push_back(names[i].substr(names[i].find("xCU[")+4));
      indices.push_back(i);
    }
  }

  assert(!x->IsHomogeneous());
  assert(indices.size());
  assert(x_space->HasIntegerMetaData("intervalID"));
  std::vector<Index> intervalIDs = x_space->GetIntegerMetaData("intervalID");

  SplitDecision retval;
  retval.parameterID = 1;
  const Number* values = x->Values();
  // boundary value is HARD coded -- 0.2
  Number reference = fabs(values[indices.at(0)]-0.2);
  retval.intervalID = intervalIDs[indices.at(0)];

  for (unsigned int i=0;i<bound_values.size();i++) {
    printf("\n%s",bound_values[i].c_str());
    bound_values[i] = bound_values[i].erase(bound_values[i].find("]"));
    printf(" turns into %s",bound_values[i].c_str());
    if (fabs(values[indices.at(i)]-0.2) < reference)
      retval.intervalID = intervalIDs[indices.at(i)];
  }

  return retval;
}

SmartPtr<DenseVector> LinearizeKKT::applyMINRESOnInterval(SmartPtr<IpoptApplication> app,const Index& interval, const Index& parameter)
{
  // printf("\n");
  // app->IpoptDataObject()->curr()->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"curr()");
  // initialize neccessary quantities
  const Index n_parameters = intervals_.getParameterCount();
  std::vector<Index> intervalIDs = intervals_.getIntervalIDVec();
  SmartPtr<DenseVectorSpace> z_space = new DenseVectorSpace(rhs_dim_);
  SmartPtr<DenseVector> z_postsplit = z_space->MakeNewDenseVector();
  SmartPtr<const DenseVector> deltap;

  SmartPtr<DenseVectorSpace> res_vspace = new DenseVectorSpace(rhs_dim_);
  SmartPtr<DenseVector> res_vsense = res_vspace->MakeNewDenseVector();

  SmartPtr<const IteratesVectorSpace> top_space;
  SmartPtr<IteratesVector> rhs_top;
  SmartPtr<DenseVector> rhs_x,rhs_s,rhs_y_c,rhs_y_d,rhs_z_L,rhs_z_U,rhs_v_L,rhs_v_U;
  SmartPtr<ShiftVector> rhs;
  SmartPtr<ShiftVector> lhs;
  SmartPtr<ShiftVector> z_cond;
  SmartPtr<ShiftVector> x0;
  Number tolo;
  Index n_m,n_rst;

  // init sensitivity matrix and -space
  SmartPtr<MultiVectorMatrixSpace> sense_space = new MultiVectorMatrixSpace(n_parameters*2,*z_space);
  SmartPtr<MultiVectorMatrix> sense = sense_space->MakeNewMultiVectorMatrix();
  SmartPtr<MultiVectorMatrix> res_sense = sense_space->MakeNewMultiVectorMatrix();

  // start algorithm on each existing intervalID
  Index s_count = 0;
  for (unsigned int i=0;i<intervalIDs.size();i++) {
    if (intervalIDs[i]==interval) {

      // split interval indices into to be shifted interval entries and remainder
      if (x_intervalIDs_.size()) {
	if (!splitIntervalIndices(x_intervalIDs_,shift_x_indices_,interval))
	  printf("\nLinearizeKKT::applyMINESOnInterval(): ERROR: unable to split x_intervalIDs_");
      } else
	printf("\nLinearizeKKT::applyMINESOnInterval(): Warning: x_intervalIDs.size() is 0.");
      if (c_intervalIDs_.size()) {
	if (!splitIntervalIndices(c_intervalIDs_,shift_c_indices_,interval))
	  printf("\nLinearizeKKT::applyMINESOnInterval(): ERROR: unable to split c_intervalIDs_");
      } else
	printf("\nLinearizeKKT::applyMINESOnInterval(): Warning: c_intervalIDs.size() is 0.");
      if (d_intervalIDs_.size()) {
	if (!splitIntervalIndices(d_intervalIDs_,shift_d_indices_,interval))
	  printf("\nLinearizeKKT::applyMINESOnInterval(): ERROR: unable to split d_intervalIDs_");
      } else
	printf("\nLinearizeKKT::applyMINESOnInterval(): Warning: d_intervalIDs.size() is 0.");

      // get different rhs entries for this column
      x_i_ = extractColumn(Wp_,i);
      y_c_i_ = extractColumn(Ap_,i);
      y_d_i_ = extractColumn(Bp_,i);

      // assign rhs subvectors
      u_i_ = shrink(x_i_,u_indices_);

      // create ShiftVector for rhs
      top_space = dynamic_cast<const IteratesVectorSpace*>(GetRawPtr(app->IpoptDataObject()->curr()->OwnerSpace()));
      rhs_top = top_space->MakeNewIteratesVector();
      rhs_top->Set(0.0);
      // x_i_->Set(5.0);
      // y_c_i_->Set(6.0);
      // y_d_i_->Set(7.0);
      rhs_top->Set_x_NonConst(*x_i_);
      rhs_top->Set_y_c_NonConst(*y_c_i_);
      rhs_top->Set_y_d_NonConst(*y_d_i_);

      SmartPtr<DenseVectorSpace> sv_bot_space = new DenseVectorSpace(int(shift_x_indices_.size()));
      rhs_x = sv_bot_space->MakeNewDenseVector();
      //      rhs_x->Set(8.0);
      rhs_x = shrink(x_i_,shift_x_indices_);

      sv_bot_space = new DenseVectorSpace(int(shift_c_indices_.size()));
      rhs_y_c = sv_bot_space->MakeNewDenseVector();
      //rhs_y_c->Set(9.0);
      rhs_y_c = shrink(y_c_i_,shift_c_indices_);

      sv_bot_space = new DenseVectorSpace(int(shift_d_indices_.size()));
      rhs_y_d = sv_bot_space->MakeNewDenseVector();
      //rhs_y_d->Set(10.0);
      rhs_y_d = shrink(y_d_i_,shift_d_indices_);

      // fill rest with 0s
      sv_bot_space = new DenseVectorSpace(int(shift_d_indices_.size()));
      rhs_s = sv_bot_space->MakeNewDenseVector();
      rhs_s->Set(0.0);

      sv_bot_space = new DenseVectorSpace(int( (Px_L_->NCols() -n_ul_) /n_i_));
      rhs_z_L = sv_bot_space->MakeNewDenseVector();
      rhs_z_L->Set(0.0);

      sv_bot_space = new DenseVectorSpace(int( (Px_U_->NCols() -n_uu_) /n_i_));
      rhs_z_U = sv_bot_space->MakeNewDenseVector();
      rhs_z_U->Set(0.0);

      sv_bot_space = new DenseVectorSpace(int(Pd_L_->NCols()/n_i_));
      rhs_v_L = sv_bot_space->MakeNewDenseVector();
      rhs_v_L->Set(0.0);

      sv_bot_space = new DenseVectorSpace(int(Pd_U_->NCols()/n_i_));
      rhs_v_U = sv_bot_space->MakeNewDenseVector();
      rhs_v_U->Set(0.0);

      rhs = new ShiftVector(rhs_top,rhs_x,rhs_s,rhs_y_c,rhs_y_d,rhs_z_L,rhs_z_U,rhs_v_L,rhs_v_U);
      // printf("\n");
      // rhs->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"rhs=tr");


      // lhs = new ShiftVector(*rhs);
      // ////////////////////////////////////////////////////////////////
      // rhs_top = top_space->MakeNewIteratesVector();
      // rhs_top->Set(0.0);
      // x_i_ = x_i_->MakeNewDenseVector();
      // x_i_->Set(31.0);
      // y_c_i_ = y_c_i_->MakeNewDenseVector();
      // y_c_i_->Set(37.0);
      // y_d_i_ = y_d_i_->MakeNewDenseVector();
      // y_d_i_->Set(41.0);
      // rhs_top->Set_x_NonConst(*x_i_);
      // rhs_top->Set_y_c_NonConst(*y_c_i_);
      // rhs_top->Set_y_d_NonConst(*y_d_i_);

      // rhs_x = x_sp->MakeNewDenseVector();
      // rhs_x->Set(43.0);
      // //      rhs_x = shrink(x_i_,shift_x_indices_);

      // rhs_y_c = y_c_sp->MakeNewDenseVector();
      // rhs_y_c->Set(47.0);
      // //      rhs_y_c = shrink(y_c_i_,shift_c_indices_);

      // rhs_y_d = y_d_sp->MakeNewDenseVector();
      // rhs_y_d->Set(53.0);
      // //      rhs_y_d = shrink(y_d_i_,shift_d_indices_);

      // lhs = new ShiftVector(rhs_top,rhs_x,rhs_y_c,rhs_y_d);

      // ////////////////////////////////////////////////////////////////
      // create x0
      x0 = new ShiftVector(*rhs);
      x0->Set(0.0);

      // assign options and compute solution of the MINRES-Step
      assignGMRESOptions(app->Options(),tolo,n_m,n_rst);
      //      testMINRES();
      //      z_cond = new ShiftVector (*x0);
      //      z_cond = computeMINRES(rhs,x0,tolo,n_m);
      rhs->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"rhs pre");
      z_cond = computeKSymMultVector(rhs);
      rhs->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"rhs pst");
      z_cond->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"KMrhs ");

      // convert symmetricized solution from MINRES back to unsymmetric equivalent
      convertBackToUnsymmetric(z_cond);

/*      printf("\n");




      // printf("\n");


      z_cond->Set(1.0);
      printf("\n");
      z_cond->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"1s");
      lhs->Set(1.0);
      printf("\n");
      lhs->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"1s");

      //      printf("\nNrm2() of z_cond is %e",z_cond->Nrm2());
      printf("\nlhs Dot x0 is %e. should be 0",lhs->Dot(*x0));
      //      printf("\nlhs Dot z_cond is %e. should be %e.",lhs->Dot(*z_cond),lhs->Nrm2()*lhs->Nrm2());

      //computeKSymMultVector();


      printf("\n");
*/
      // printf("\n");
      // rhs->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"rhs");
      // printf("\n");
      // lhs->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"lhs");
      // x0 = computeAMultVector(z_cond);
      // x0->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"prhs");

      // rhs->AddOneVector(-1.0,*x0,1.0);
      // printf("\n");

      // printf("\n");

      // x0 = computeKIMultVector(lhs);
      // x0->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"Alhs");
      // z_cond = computeAMultVector(x0);
      // z_cond->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"plhs");

      // lhs->AddOneVector(-1.0,*z_cond,1.0);
      // printf("\n");
      // lhs->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"other 0s?");

      // convert the ShiftVector into a DenseVector
      z_postsplit = z_cond->getDVector();

      // add the resulting sense vector to matrix ....
      sense->SetVector(s_count,*z_postsplit);
      s_count++;
    }
  }
  // get the parameter perturbation for this split
  deltap = computeDeltaP(intervals_,interval,parameter);
  // calculate absolute sensitivity values
  sense->MultVector(1.0,*deltap,0.0,*res_vsense);

  return res_vsense;
}

std::vector<Index> LinearizeKKT::collectIntMetaData() const
{
  //  printf("\nLinearizeKKT::collectIntMetaData(): rhs_dim_ is: %d",rhs_dim_);
  std::vector<Index> retval(rhs_dim_);
  Index abs_pos = 0;
  // x top
  //printf("\nLinearizeKKT::collectIntMetaData(): abs_pos is: %d    x_intervalIDs_.size() is: %d  ",abs_pos,int(x_intervalIDs_.size()));
  for (unsigned int m=0;m<x_intervalIDs_.size();m++) {
    retval[abs_pos] = x_intervalIDs_[m];
    abs_pos++;
  }
  // s top
//  printf("\nLinearizeKKT::collectIntMetaData(): abs_pos is: %d    d_intervalIDs_.size() is: %d  ",abs_pos,int(d_intervalIDs_.size()));
  for (unsigned int m=0;m<d_intervalIDs_.size();m++) {
    retval[abs_pos] = d_intervalIDs_[m];
    abs_pos++;
  }
  // y_c top
//  printf("\nLinearizeKKT::collectIntMetaData(): abs_pos is: %d    c_intervalIDs_.size() is: %d  ",abs_pos,int(c_intervalIDs_.size()));
  for (unsigned int m=0;m<c_intervalIDs_.size();m++) {
    retval[abs_pos] = c_intervalIDs_[m];
    abs_pos++;
  }
  // y_d top
//  printf("\nLinearizeKKT::collectIntMetaData(): abs_pos is: %d    d_intervalIDs_.size() is: %d  ",abs_pos,int(d_intervalIDs_.size()));
  for (unsigned int m=0;m<d_intervalIDs_.size();m++) {
    retval[abs_pos] = d_intervalIDs_[m];
    abs_pos++;
  }
  // z_L top
  SmartPtr<const ExpansionMatrix> exp = dynamic_cast<const ExpansionMatrix*>(GetRawPtr(Px_L_));
  const Index* epi_zl = exp->ExpandedPosIndices();
//  printf("\nLinearizeKKT::collectIntMetaData(): abs_pos is: %d    exp->NCols() is: %d  ",abs_pos,exp->NCols());
  for (int m=0;m<exp->NCols();m++) {
    retval[abs_pos] = x_intervalIDs_[*(epi_zl+m)];
    abs_pos++;
  }
  // z_U top
  exp = dynamic_cast<const ExpansionMatrix*>(GetRawPtr(Px_U_));
  const Index* epi_zu = exp->ExpandedPosIndices();
//  printf("\nLinearizeKKT::collectIntMetaData(): abs_pos is: %d    exp->NCols() is: %d  ",abs_pos,exp->NCols());
  for (int m=0;m<exp->NCols();m++) {
    retval[abs_pos] = x_intervalIDs_[*(epi_zu+m)];
    abs_pos++;
  }
  // v_L top
  exp = dynamic_cast<const ExpansionMatrix*>(GetRawPtr(Pd_L_));
  const Index* epi_vl = exp->ExpandedPosIndices();
//  printf("\nLinearizeKKT::collectIntMetaData(): abs_pos is: %d    exp->NCols() is: %d  ",abs_pos,exp->NCols());
  for (int m=0;m<exp->NCols();m++) {
    retval[abs_pos] = d_intervalIDs_[*(epi_vl+m)];
    abs_pos++;
  }
  // v_U top
  exp = dynamic_cast<const ExpansionMatrix*>(GetRawPtr(Pd_U_));
  const Index* epi_vu = exp->ExpandedPosIndices();
//  printf("\nLinearizeKKT::collectIntMetaData(): abs_pos is: %d    exp->NCols() is: %d  ",abs_pos,exp->NCols());
  for (int m=0;m<exp->NCols();m++) {
    retval[abs_pos] = d_intervalIDs_[*(epi_vu+m)];
    abs_pos++;
  }

  // x shift
  for (unsigned int m=0;m<shift_x_indices_.size();m++) {
    retval[abs_pos] = intervals_.size();
    abs_pos++;
  }
  // s shift
  for (unsigned int m=0;m<shift_d_indices_.size();m++) {
    retval[abs_pos] = intervals_.size();
    abs_pos++;
  }
  // y_c shift
  for (unsigned int m=0;m<shift_c_indices_.size();m++) {
    retval[abs_pos] = intervals_.size();
    abs_pos++;
  }
  // y_d shift
  for (unsigned int m=0;m<shift_d_indices_.size();m++) {
    retval[abs_pos] = intervals_.size();
    abs_pos++;
  }
  // z_L shift
  Index dim = int((z_L_->Dim() - (z_L_->Dim() % n_i_)) / n_i_);
  for (int m=0;m<dim;m++) {
    retval[abs_pos] = intervals_.size();
    abs_pos++;
  }
  // z_U shift
  dim = int((z_U_->Dim() - (z_U_->Dim() % n_i_)) / n_i_);
  for (int m=0;m<dim;m++) {
    retval[abs_pos] = intervals_.size();
    abs_pos++;
  }
  // v_L shift
  dim = int((v_L_->Dim() - (v_L_->Dim() % n_i_)) / n_i_);
  for (int m=0;m<dim;m++) {
    retval[abs_pos] = intervals_.size();
    abs_pos++;
  }
  // v_U shift
  dim = int((v_U_->Dim() - (v_U_->Dim() % n_i_)) / n_i_);
  for (int m=0;m<dim;m++) {
    retval[abs_pos] = intervals_.size();
    abs_pos++;
  }

  return retval;
}

SmartPtr<ShiftVector> LinearizeKKT::computeMINRES(SmartPtr<ShiftVector>b,SmartPtr<ShiftVector>x0,const Number& tol, const Index& n_max, const Index& n_rest) const
{
  std::vector<SmartPtr<ShiftVector> > r(2);
  std::vector<SmartPtr<ShiftVector> > s(3);
  std::vector<SmartPtr<ShiftVector> > p(3);
  Number alpha;
  std::vector<Number> beta(2);
  std::vector<SmartPtr<ShiftVector> > x(2);
  Index j_thresh = 1;
  SmartPtr<ShiftVector> retval;

  x[1] = new ShiftVector(*x0);
  r[1] = new ShiftVector(*b);
  r[1]->AddOneVector(-1.0,*computeKSymMultVector(x[1]),1.0);
  p[1] = new ShiftVector(*r[1]);
  p[0] = new ShiftVector(*p[1]);
  s[1] = computeKSymMultVector(p[1]);
  s[0] = new ShiftVector(*s[1]);
  // printf("\n");
  // s[0]->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"s[0]");

  std::string fname = "minres.dat";
  std::ofstream minres;
  char buffer[63];
  minres.open(fname.c_str());
  minres << "#.dat file automatically generated by AMPL intervallization routine\n#Ben Waldecker Sep 2012\n";

  for (int i=1;i<n_max;i++) {
    // recursive value shift
    p[2] = new ShiftVector(*p[1]);
    p[1] = new ShiftVector(*p[0]);
    s[2] = new ShiftVector(*s[1]);
    s[1] = new ShiftVector(*s[0]);
    r[0] = new ShiftVector(*r[1]);
    x[0] = new ShiftVector(*x[1]);

    if (s[1] ->Dot(* s[1]) == 0) {
    printf("\nLinearizeKKT:computeMINRES(): Assertion due to 's' being a douche: Nrm2 = %e  or %f",s[1]->Nrm2(),s[1]->Nrm2());
      s[1]->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"s_assert");
    }

    assert(s[1] ->Dot(* s[1]));

    // calculate steplength
    alpha = (r[0] ->Dot(* s[1])) / (s[1] ->Dot(* s[1]));
    // update solution
    x[1] = new ShiftVector(*x[0]);
    x[1]->AddOneVector(alpha,*p[1],1.0);

    // output current residual
    SmartPtr<ShiftVector> tmp = computeKSymMultVector(x[1]);
    tmp->AddOneVector(-1.0,*b,1.0);
    sprintf(buffer,"\nresidual[%d] = %e     (relative: %e)",i,tmp->Nrm2(), tmp->Nrm2()/b->Nrm2());
    minres << buffer;

    //update residual
    r[1] = new ShiftVector(*r[0]);
    r[1]->AddOneVector(-1.0*alpha,*s[1],1.0);

    if (r[1]->Nrm2()<tol) {

      // target accuracy achieved - compute and return solution
      printf("\n target accuracy achieved. r[1]->Nrm2() is: %e, which is less than tol: %e\n",r[1]->Nrm2(),tol);
      x[1]->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"x_sol");

      retval = new ShiftVector(*x[1]);
      return retval;

    } else {
      // solution needs to improve - update next stepdirection
      p[0] = new ShiftVector(*s[1]);

      // update s
      s[0] = computeKSymMultVector(s[1]);

      if (i>=2)
	j_thresh=2;

      // modify updated  p and s
    assert(s[2] ->Dot(* s[2]));
      for (int j=0;j<j_thresh;j++) {
	beta[j] = s[j+1] ->Dot(* s[0] ) / s[j+1] ->Dot(* s[j+1]);
	p[0]->AddOneVector(-1.0*beta[j],*p[j+1],1.0);
	s[0]->AddOneVector(-1.0*beta[j],*s[j+1],1.0);
      }
    }
  }

  minres << "\n\n#end of file";
  minres.close();

  return new ShiftVector(*x[1]);
}

void LinearizeKKT::convertBackToUnsymmetric(SmartPtr<ShiftVector> target) const
{
  // hard coded - needs to be read from options
  // const Number mu_target = readMuFromOptions(options_);
  // ODER in GMRES Options auch assignen
  const Number mu_target = 1.0e-6;

  // top - z
  SmartPtr<DenseVector> top_zl = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  SmartPtr<DenseVector> top_x = dynamic_cast<DenseVector*>(target->top()->x()->MakeNewCopy());
  SmartPtr<DenseVector> z_B_e = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  top_zl->Set(mu_target);
  top_zl->ElementWiseDivide(*Sl_x_L_);
  top_x->ElementWiseMultiply(*z_L_DIV_Sl_);
  Px_L_->MultVector(1.0,*z_L_,0.0,*z_B_e);
  z_B_e->AddTwoVectors(-1.0,*top_x,1.0,*top_zl,-1.0);
  top_zl = dynamic_cast<DenseVector*>(target->top()->z_L()->MakeNew());
  Px_L_->TransMultVector(1.0,*z_B_e,0.0,*top_zl);

  SmartPtr<DenseVector> top_zu = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  top_x = dynamic_cast<DenseVector*>(target->top()->x()->MakeNewCopy());
  z_B_e = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  top_zu->Set(mu_target);
  top_zu->ElementWiseDivide(*Sl_x_U_);
  top_x->ElementWiseMultiply(*z_U_DIV_Sl_);
  Px_U_->MultVector(1.0,*z_U_,0.0,*z_B_e);
  z_B_e->AddTwoVectors(-1.0,*top_x,1.0,*top_zu,-1.0);
  top_zu = dynamic_cast<DenseVector*>(target->top()->z_U()->MakeNew());
  Px_U_->TransMultVector(1.0,*z_B_e,0.0,*top_zu);

  // top - v
  SmartPtr<DenseVector> top_vl = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  SmartPtr<DenseVector> top_s = dynamic_cast<DenseVector*>(target->top()->s()->MakeNewCopy());
  SmartPtr<DenseVector> v_B_e = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  top_vl->Set(mu_target);
  top_vl->ElementWiseDivide(*Sl_s_L_);
  top_s->ElementWiseMultiply(*v_L_DIV_Sl_);
  Pd_L_->MultVector(1.0,*v_L_,0.0,*v_B_e);
  v_B_e->AddTwoVectors(-1.0,*top_s,1.0,*top_vl,-1.0);
  top_vl = dynamic_cast<DenseVector*>(target->top()->v_L()->MakeNew());
  Pd_L_->TransMultVector(1.0,*v_B_e,0.0,*top_vl);

  SmartPtr<DenseVector> top_vu = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  top_s = dynamic_cast<DenseVector*>(target->top()->s()->MakeNewCopy());
  v_B_e = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  top_vu->Set(mu_target);
  top_vu->ElementWiseDivide(*Sl_s_U_);
  top_s->ElementWiseMultiply(*v_U_DIV_Sl_);
  Pd_U_->MultVector(1.0,*v_U_,0.0,*v_B_e);
  v_B_e->AddTwoVectors(-1.0,*top_s,1.0,*top_vu,-1.0);
  top_vu = dynamic_cast<DenseVector*>(target->top()->v_U()->MakeNew());
  Pd_U_->TransMultVector(1.0,*v_B_e,0.0,*top_vu);

  // bot - z
  SmartPtr<DenseVector> bot_zl = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  SmartPtr<DenseVector> bot_x;
  z_B_e = dynamic_cast<DenseVector*>(target->x()->MakeNewCopy());
  // printf("\n");
  // z_B_e->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"zbe pre");
  bot_x = expand(z_B_e,shift_x_indices_,bot_zl->Dim());
  bot_zl->Set(mu_target);
  bot_zl->ElementWiseDivide(*Sl_x_L_);
  bot_x->ElementWiseMultiply(*z_L_DIV_Sl_);
  // printf("\n");
  // z_B_e->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"zbe pst");
  z_B_e = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  Px_L_->MultVector(1.0,*z_L_,0.0,*z_B_e);
  z_B_e->AddTwoVectors(-1.0,*bot_x,1.0,*bot_zl,-1.0);
  bot_zl = shrink(z_B_e,shift_x_indices_);
  // printf("\n");
  // bot_zl->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"bzl pre");

  SmartPtr<DenseVector> bot_zu = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  z_B_e = dynamic_cast<DenseVector*>(target->x()->MakeNewCopy());
  // printf("\n");
  // bot_zl->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"bzl pst");
  bot_x = expand(z_B_e,shift_x_indices_,bot_zu->Dim());
  bot_zu->Set(mu_target);
  bot_zu->ElementWiseDivide(*Sl_x_U_);
  bot_x->ElementWiseMultiply(*z_U_DIV_Sl_);
  z_B_e = dynamic_cast<DenseVector*>(target->top()->x()->MakeNew());
  Px_U_->MultVector(1.0,*z_U_,0.0,*z_B_e);
  z_B_e->AddTwoVectors(-1.0,*bot_x,1.0,*bot_zu,-1.0);
  bot_zu = shrink(z_B_e,shift_x_indices_);

  // bot - v
  SmartPtr<DenseVector> bot_vl = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  SmartPtr<DenseVector> bot_s;
  v_B_e = dynamic_cast<DenseVector*>(target->s()->MakeNewCopy());
  bot_s = expand(v_B_e,shift_d_indices_,bot_vl->Dim());
  bot_vl->Set(mu_target);
  bot_vl->ElementWiseDivide(*Sl_s_L_);
  bot_s->ElementWiseMultiply(*v_L_DIV_Sl_);
  v_B_e = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  Pd_L_->MultVector(1.0,*v_L_,0.0,*v_B_e);
  v_B_e->AddTwoVectors(-1.0,*bot_s,1.0,*bot_vl,-1.0);
  bot_vl = shrink(v_B_e,shift_d_indices_);

  SmartPtr<DenseVector> bot_vu = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  v_B_e = dynamic_cast<DenseVector*>(target->s()->MakeNewCopy());
  bot_s = expand(v_B_e,shift_d_indices_,bot_vu->Dim());
  bot_vu->Set(mu_target);
  bot_vu->ElementWiseDivide(*Sl_s_U_);
  bot_s->ElementWiseMultiply(*v_U_DIV_Sl_);
  v_B_e = dynamic_cast<DenseVector*>(target->top()->s()->MakeNew());
  Pd_U_->MultVector(1.0,*v_U_,0.0,*v_B_e);
  v_B_e->AddTwoVectors(-1.0,*bot_s,1.0,*bot_vu,-1.0);
  bot_vu = target->v_U()->MakeNewDenseVector();
  bot_vu = shrink(v_B_e,shift_d_indices_);

  // insert results into target
  target->Set_z_L_top(*top_zl);
  target->Set_z_U_top(*top_zu);
  target->Set_v_L_top(*top_vl);
  target->Set_v_U_top(*top_vu);
  target->Set_z_L(*bot_zl);
  target->Set_z_U(*bot_zu);
  target->Set_v_L(*bot_vl);
  target->Set_v_U(*bot_vu);

  // printf("\n");
  // target->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"cBTU");
}

void LinearizeKKT::testGMRES()
{
  // initialize neccessary vars
  SmartPtr<DenseVectorSpace> space = new DenseVectorSpace(80);
  SmartPtr<DenseVector> b = space->MakeNewDenseVector();
  SmartPtr<MultiVectorMatrixSpace> m_space = new MultiVectorMatrixSpace(80,*space);
  SmartPtr<MultiVectorMatrix> A = m_space->MakeNewMultiVectorMatrix();

  // setup A
  Index abs_pos = 0;
  Number* vals = new Number[80];
  for (int i=0;i<10;i++) {
    SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
    for (int j=0;j<80;j++) {
      vals[j] = j+i*i;
    }
    moo->SetValues(vals);
    A->SetVector(abs_pos,*moo);
    abs_pos++;
  }
  for (int i=0;i<10;i++) {
    SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
    for (int j=0;j<80;j++) {
      vals[j] = j*i+i;
    }
    moo->SetValues(vals);
    A->SetVector(abs_pos,*moo);
    abs_pos++;
  }
  for (int i=0;i<10;i++) {
    SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
    for (int j=0;j<40;j++) {
      vals[2*j] = j;
      vals[2*j+1] = 0;
    }
    moo->SetValues(vals);
    A->SetVector(abs_pos,*moo);
    abs_pos++;
  }  for (int i=0;i<10;i++) {
    SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
    for (int j=0;j<40;j++) {
      vals[2*j] = 0 ;
      vals[2*j+1] = j;
    }
    moo->SetValues(vals);
    A->SetVector(abs_pos,*moo);
    abs_pos++;
  }  for (int i=0;i<10;i++) {
    SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
    for (int j=0;j<40;j++) {
      vals[2*j] = 1;
      vals[2*j+1] = 0;
    }
    moo->SetValues(vals);
    A->SetVector(abs_pos,*moo);
    abs_pos++;
  }  for (int i=0;i<10;i++) {
    SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
    for (int j=0;j<40;j++) {
      vals[2*j] = 0;
      vals[2*j+1] = -2;
    }
    moo->SetValues(vals);
    A->SetVector(abs_pos,*moo);
    abs_pos++;
  }  for (int i=0;i<10;i++) {
    SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
    for (int j=0;j<80;j++) {
      if (j<40)
	vals[j] = 1;
      else
	vals[j] = 0;
    }
    moo->SetValues(vals);
    A->SetVector(abs_pos,*moo);
    abs_pos++;
  }  for (int i=0;i<10;i++) {
    SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
    for (int j=0;j<80;j++) {
      if (j<40)
	vals[j] = 0;
      else
	vals[j] = 1;
    }
    moo->SetValues(vals);
    A->SetVector(abs_pos,*moo);
    abs_pos++;
  }

  // setup b
  for (int j=0;j<80;j++) {
        vals[j] = 80-j;
  }
  b->SetValues(vals);
  A->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"testgmres: A");
  Index n_max = 81;
  SmartPtr<DenseVector> r = b;
  SmartPtr<DenseVector> x0 = space->MakeNewDenseVector();
  x0->Set(0);
  SmartPtr<DenseVector> x = dynamic_cast<DenseVector*>(x0->MakeNewCopy());

  std::vector<SmartPtr<DenseVector> > v(n_max+1);
  std::vector<SmartPtr<DenseVector> > w(n_max+1);
  Number tol = 1e-9;
  Number beta;
  Index g_cnt = 0;
  std::vector<Number> h((n_max+1)*(n_max+1));
  std::vector<Number> s(n_max+1);
  std::vector<Number> c(n_max+1);
  std::vector<Number> gamma(n_max+1);
  std::vector<Number> y(n_max);

  SmartPtr<DenseVector> tmp;

  // for residual output
  std::string fname = "tresiduals.dat";
  std::ofstream residuals;
  residuals.open(fname.c_str());
  char buffer[63];
  residuals << "#.dat file automatically generated by AMPL intervallization routine\n#Ben Waldecker Sep 2012\n";

  // preloop inital var values
  //  printf("\nnorm of r is %f",r->Nrm2());
  v[0] = dynamic_cast<DenseVector*>(r->MakeNewCopy());
  gamma[0] = r->Nrm2();
  assert(gamma[0]);
  v[0]->Scal(1/gamma[0]);
  printf("\ngamma[0] is %f. Nrm of v[0] is %e",r->Nrm2(),v[0]->Nrm2());

  // outer loop
  for (Index j=0; j<n_max;j++) {
  //  printf("\ncomputeGMRES(): outer loop iteration started. j = %d",j);

    //w[j] = computeKPMultVector(v[j]);
    w[j] = space->MakeNewDenseVector();
    A->MultVector(1.0,*v[j],0.0,*w[j]);

    // update h[i,j] entries
    for (Index i=0; i<=j;i++) {
      h[i+j*n_max] = v[i]->Dot(*w[j]);
    }

    // calculate w[j]
    for (Index i=0;i<=j;i++) {
      w[j]->AddOneVector(-1.0*h[i+j*n_max],*v[i],1.0);
    }
    // rotate it
    for (Index i=0;i<j;i++) {
      double hij = h[i+j*n_max];
      double hij1 = h[i+1+j*n_max];
      h[i+j*n_max] = c[i+1]*hij + s[i+1]*hij1;
      h[i+1+j*n_max] = -s[i+1]*hij + c[i+1]*hij1;
    }
    // expand h, first step
    h[j+1+j*n_max] = w[j]->Nrm2();

    // expand h, second step - eigenvector style

    beta = sqrt(h[j+j*n_max]*h[j+j*n_max]+h[j+1+j*n_max]*h[j+1+j*n_max]);
    s[j+1] = h[j+1+j*n_max]/beta;
    c[j+1] = h[j+j*n_max]/beta;
    h[j+j*n_max] = beta;

    // expand and update gammas
    gamma[j+1] = -s[j+1]*gamma[j];
    gamma[j] = c[j+1]*gamma[j];
    g_cnt++;

    // check quality of current solution
    if ( ( (fabs(gamma[j+1]) > tol) && (!n_max) ) || ( (fabs(gamma[j+1]) > tol) && n_max && (j!=n_max-1) ))  {
      // insufficient - update v[j+1]
      assert(h[j+1+j*n_max]);
      v[j+1] = dynamic_cast<DenseVector*>(w[j]->MakeNewCopy());
      v[j+1]->Scal(1/h[j+1+j*n_max]);
   } else {
     // set up multipliers y[i]
     for (Index i=j;i>=0;i--) {
       assert(h[i+i*n_max]);

       // set up sum[k] h[i+k*n_max] y[k]
       beta = 0;
       for (Index k=i+1;k<=j;k++) {
	 beta += h[i+k*n_max]*y[k];
       }
       y[i]= (1/h[i+i*n_max])*(gamma[i]-beta);
     }

     // compute solution
     for (Index i=0; i<=j;i++) {
       x->AddOneVector(y[i],*v[i],1.0);
     }
     break;
   }

    /* set up current solution to calculate absolute residuum and write it to file
       in order to check for convergence */
      for (Index i=j;i>=0;i--) {
	assert(h[i+i*n_max]);
	beta = 0;
	for (Index k=i+1;k<=j;k++) {
	  beta += h[i+k*n_max]*y[k];
	}
	y[i]= 1/h[i+i*n_max]*(gamma[i]-beta);
      }
      for (Index i=0; i<=j;i++) {
	x->AddOneVector(y[i],*v[i],1.0);
      }
      SmartPtr<DenseVector> fooo = space->MakeNewDenseVector();
      A->MultVector(1.0,*x,0.0,*fooo);
      fooo->AddOneVector(-1.0,*b,1.0);
    Number residual = fooo->Nrm2();

    //write residual into residuals.dat file for python access
    sprintf(buffer,"\nresidual[%d] = %e, relative residual: %e",j,residual, residual/b->Nrm2());
    residuals << buffer;
    x = dynamic_cast<DenseVector*>(x0->MakeNewCopy());
  }
  residuals << "\n\n#end of file";
  residuals.close();

  //write gamma[j] into gamma.dat file for python access
  fname = "tgammas.dat";
  std::ofstream gammas;
  gammas.open(fname.c_str());
  gammas << "#.dat file automatically generated by AMPL intervallization routine\n#Ben Waldecker Sep 2012\n";
  for (int i=0;i<g_cnt;i++) {
    sprintf(buffer,"\ngamma[%d] = %e",i,gamma[i]);
    gammas << buffer;
  }
  gammas << "\n\n#end of file";
  gammas.close();

}

void LinearizeKKT::testMINRES()
{
  /**/
  // initialize neccessary vars
  SmartPtr<DenseVectorSpace> space = new DenseVectorSpace(8);
  SmartPtr<DenseVector> b = space->MakeNewDenseVector();
  SmartPtr<MultiVectorMatrixSpace> m_space = new MultiVectorMatrixSpace(8,*space);
  SmartPtr<MultiVectorMatrix> A = m_space->MakeNewMultiVectorMatrix();

  // setup A
  Index abs_pos = 0;
  Number* vals0 = new Number[8];
  Number* vals1 = new Number[8];
  Number* vals2 = new Number[8];
  Number* vals3 = new Number[8];
  Number* vals4 = new Number[8];
  Number* vals5 = new Number[8];
  Number* vals6 = new Number[8];
  Number* vals7 = new Number[8];

  vals0[0] = 1.0; vals1[0] = 0.0; vals2[0] = 0.0; vals3[0] = 3.0; vals4[0] = 0.0; vals5[0] = 0.0; vals6[0] = 0.0; vals7[0] = 0.0;
  vals0[1] = 0.0; vals1[1] = 1.0; vals2[1] = 0.0; vals3[1] = 0.0; vals4[1] = 4.0; vals5[1] = 0.0; vals6[1] = 0.0; vals7[1] = 0.0;
  vals0[2] = 0.0; vals1[2] = 0.0; vals2[2] =-2.0; vals3[2] = 0.0; vals4[2] = 0.0; vals5[2] = 1.0; vals6[2] = 0.0; vals7[2] = 0.0;
  vals0[3] = 3.0; vals1[3] = 0.0; vals2[3] = 0.0; vals3[3] =-1.0; vals4[3] = 0.0; vals5[3] = 0.0; vals6[3] = 0.0; vals7[3] = 0.0;
  vals0[4] = 0.0; vals1[4] = 4.0; vals2[4] = 0.0; vals3[4] = 0.0; vals4[4] =-1.0; vals5[4] = 0.0; vals6[4] = 0.0; vals7[4] = 0.0;
  vals0[5] = 0.0; vals1[5] = 0.0; vals2[5] = 1.0; vals3[5] = 0.0; vals4[5] = 0.0; vals5[5] =-2.0; vals6[5] = 0.0; vals7[5] = 0.0;
  vals0[6] = 0.0; vals1[6] = 0.0; vals2[6] = 0.0; vals3[6] = 0.0; vals4[6] = 0.0; vals5[6] = 0.0; vals6[6] = 4.0; vals7[6] = 0.0;
  vals0[7] = 0.0; vals1[7] = 0.0; vals2[7] = 0.0; vals3[7] = 0.0; vals4[7] = 0.0; vals5[7] = 0.0; vals6[7] = 0.0; vals7[7] = 3.0;

  SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
  moo->SetValues(vals0);
  A->SetVector(abs_pos,*moo);
  abs_pos++;
  moo = space->MakeNewDenseVector();
  moo->SetValues(vals1);
  A->SetVector(abs_pos,*moo);
  abs_pos++;
  moo = space->MakeNewDenseVector();
  moo->SetValues(vals2);
  A->SetVector(abs_pos,*moo);
  abs_pos++;
  moo = space->MakeNewDenseVector();
  moo->SetValues(vals3);
  A->SetVector(abs_pos,*moo);
  abs_pos++;
  moo = space->MakeNewDenseVector();
  moo->SetValues(vals4);
  A->SetVector(abs_pos,*moo);
  abs_pos++;
  moo = space->MakeNewDenseVector();
  moo->SetValues(vals5);
  A->SetVector(abs_pos,*moo);
  abs_pos++;
  moo = space->MakeNewDenseVector();
  moo->SetValues(vals6);
  A->SetVector(abs_pos,*moo);
  abs_pos++;
  moo = space->MakeNewDenseVector();
  moo->SetValues(vals7);
  A->SetVector(abs_pos,*moo);

  // setup b
  for (int j=0;j<8;j++) {
    *(vals1+j) = 10-j;
  }
   /*
  SmartPtr<DenseVectorSpace> space = new DenseVectorSpace(2);
  SmartPtr<DenseVector> b = space->MakeNewDenseVector();
  SmartPtr<MultiVectorMatrixSpace> m_space = new MultiVectorMatrixSpace(2,*space);
  SmartPtr<MultiVectorMatrix> A = m_space->MakeNewMultiVectorMatrix();

  // setup A
  Index abs_pos = 0;
  Number* vals0 = new Number[2];
  Number* vals1 = new Number[2];

  vals0[0] = 1.0; vals1[0] = 2.0;
  vals0[1] = 2.0; vals1[1] =-1.0;

  SmartPtr<DenseVector> moo = space->MakeNewDenseVector();
  moo->SetValues(vals0);
  A->SetVector(0,*moo);

  moo = space->MakeNewDenseVector();
  moo->SetValues(vals1);
  A->SetVector(1,*moo);
*/
  //setup b
  for (int j=0;j<2;j++) {
    *(vals1+j) = 10-j;
  }
  /**/
  b->SetValues(vals1);
  // printf("\n");
  // A->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"minres: A");
  // printf("\n");
  // b->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"minres: b");

  Number tol = 1.0e-10;
  Index n_max = 400;
  std::vector<SmartPtr<DenseVector> > r(2);
  std::vector<SmartPtr<DenseVector> > s(3);
  std::vector<SmartPtr<DenseVector> > p(3);
  Number alpha;
  std::vector<Number> beta(2);
  std::vector<SmartPtr<DenseVector> > x(2);
  Index j_thresh = 1;

  x[1] = space->MakeNewDenseVector();
  x[1]->Set(0.0);
  r[1] = dynamic_cast<DenseVector*>(b->MakeNewCopy());
  SmartPtr<DenseVector> tmp = space->MakeNewDenseVector();
  A->MultVector(1.0,*x[1],0.0,*tmp);
  r[1]->AddOneVector(-1.0,*tmp,1.0);
  p[0] = dynamic_cast<DenseVector*>(r[1]->MakeNewCopy());
  p[1] = p[0];
  s[0] = space->MakeNewDenseVector();
  A->MultVector(1.0,*p[0],0.0,*s[0]);
  s[1] = s[0];
  // printf("\n");
  // s[0]->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"s[0]=Ab");

  std::string fname = "tminres.dat";
  std::ofstream minres;
  char buffer[63];
  minres.open(fname.c_str());
  minres << "#.dat file automatically generated by AMPL intervallization routine\n#Ben Waldecker Sep 2012\n";

  for (int i=1;i<n_max;i++) {
    // recursive value shift
    p[2] = p[1];
    p[1] = p[0];
    s[2] = s[1];
    s[1] = s[0];
    r[0] = r[1];
    x[0] = x[1];

    // calculate steplength
    alpha = ( r[0] ->Dot(* s[1]) ) / ( s[1] ->Dot(* s[1] ));
    // update solution
    x[1] = dynamic_cast<DenseVector*>(x[0]->MakeNewCopy());
    x[1]->AddOneVector(alpha,*p[1],1.0);

    // output current residual
    tmp = space->MakeNewDenseVector();
    A->MultVector(1.0,*x[1],0.0,*tmp);
    tmp->AddOneVector(-1.0,*b,1.0);
    sprintf(buffer,"\nresidual[%d] = %e     (relative: %e)",i,tmp->Nrm2(), tmp->Nrm2()/b->Nrm2());
    minres << buffer;

    //update residual
    r[1] = dynamic_cast<DenseVector*>(r[0]->MakeNewCopy());
    r[1]->AddOneVector(-1.0*alpha,*s[1],1.0);

    if (r[1]->Nrm2()<tol) {
      // target accuracy achieved - compute and return solution
      printf("\n target accuracy achieved. r[1]->Nrm2() is: %e, which is less than tol: %e\n",r[1]->Nrm2(),tol);
      x[1]->Print(*jnl_, J_INSUPPRESSIBLE, J_DBG,"x_sol");
      break;
    } else {
      // solution needs to improve - update next stepdirection
      p[0] = dynamic_cast<DenseVector*>(s[1]->MakeNewCopy());

      // update s
      s[0] = space->MakeNewDenseVector();
      A->MultVector(1.0,*s[1],0.0,*s[0]);

      if (i>=2)
	j_thresh=2;

      // modify updated  p and s
      for (int j=0;j<j_thresh;j++) {
	beta[j] = s[j+1] ->Dot(* s[0] ) / s[j+1] ->Dot(* s[j+1]);
	p[0]->AddOneVector(-1.0*beta[j],*p[j+1],1.0);
	s[0]->AddOneVector(-1.0*beta[j],*s[j+1],1.0);
      }
    }
  }

  minres << "\n\n#end of file";
  minres.close();

}
