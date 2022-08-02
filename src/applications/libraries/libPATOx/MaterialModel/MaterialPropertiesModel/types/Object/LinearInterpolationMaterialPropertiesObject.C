#include "LinearInterpolationMaterialPropertiesObject.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LinearInterpolationMaterialPropertiesObject::LinearInterpolationMaterialPropertiesObject
(
    RectangularMatrix<scalar>& materialPropertyTable
)
  :
materialPropertyTable_(materialPropertyTable),
indexPiMin_(0),
indexPiMax_(0),
indexPiPlus1Min_(0),
indexPiPlus1Max_(0),
f_p_(0.0),
indexTi_(0.0),
indexTiPlus1_(0.0),
f_T_(0.0),
cp_(0.0),
h_(0.0),
ki_(0.0),
kj_(0.0),
kk_(0.0),
eps_(0.0),
alpha_(0.0),
T_temp_(0.0),
p_temp_(0.0),
small_(0.0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LinearInterpolationMaterialPropertiesObject::~LinearInterpolationMaterialPropertiesObject()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// indexTi: Returns the Temperature interpolation factor and updates indexTi_ and indexTiPlus1_
Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::indexTi(scalar T_, int i_min_, int i_max_)
{

  if
  (
      (T_ - materialPropertyTable_[i_min_][1] < -small_) ||
      (T_ - materialPropertyTable_[i_max_][1] > small_)
  ) {
//        Info << "Temperature out of materialPropertyTable range" << nl;
    if ((T_ - materialPropertyTable_[i_min_][1]) < - small_) {
//            Info << "T = "
//                 << T_
//                 << " < "
//                 << materialPropertyTable_[i_min_][1]
//                 << nl;
      T_ =   materialPropertyTable_[i_min_][1];
    }

    if ((T_ - materialPropertyTable_[i_max_][1]) > small_) {
//            Info << "T = "
//                 << T_
//                 << " > "
//                 << materialPropertyTable_[i_max_][1]
//                 << nl;
      T_ =  materialPropertyTable_[i_max_][1];
    }

  }

  int i_Ti(i_min_);

  while (materialPropertyTable_[i_Ti][1] < T_) {
    i_Ti++;
  }

  indexTi_ = i_Ti - 1;

  if (indexTi_<0) {
    indexTi_=0;
  }

  if (T_ == materialPropertyTable_[indexTi_][1]) {
    indexTiPlus1_ = indexTi_;
    f_T_ = 1;
    return f_T_;
  } else {
    indexTiPlus1_ = indexTi_ + 1;
    f_T_ =
        (
            T_ -
            materialPropertyTable_[indexTi_][1]
        ) /
        (
            materialPropertyTable_[indexTiPlus1_][1] -
            materialPropertyTable_[indexTi_][1]
        );
  }

  // Info << "f_T = " << f_T_ << nl;
  return f_T_;
}

// indexPi: Returns the interpolation factor f and updates IndexPiMin_, IndexPiMax_, IndexPiPlus1Min_, IndexPiPlus1Max_
Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::indexPi(scalar p_)
{
  int i_Pi(0);

  if
  (
      (p_ - materialPropertyTable_[0][0] > small_ ) ||
      (p_ - materialPropertyTable_[materialPropertyTable_.m() - 1][0] < -small_)
  ) {
//        Info << "Pressure out of materialPropertyTable range" << nl;

    if ((p_ - materialPropertyTable_[0][0]) > small_ ) {
//            Info << "p = "
//                 << p_
//                 << " > "
//                 << materialPropertyTable_[0][0]
//                 << nl;
      p_ = materialPropertyTable_[0][0];
    }

    if ((p_ - materialPropertyTable_[materialPropertyTable_.m() - 1][0]) < - small_) {
//            Info << "p = "
//                 << p_
//                 << " < "
//                 << materialPropertyTable_[materialPropertyTable_.m() - 1][0]
//                 << nl;
      p_= materialPropertyTable_[materialPropertyTable_.m() - 1][0];
    }

  }

  while (materialPropertyTable_[i_Pi][0] > p_) {
    i_Pi++;
  }
  indexPiMin_ = i_Pi;

  while
  (
      materialPropertyTable_[i_Pi][0] ==
      materialPropertyTable_[indexPiMin_][0]
  ) {
    i_Pi++;
    if (i_Pi == materialPropertyTable_.m()) {
      break;
    }
  }

  indexPiMax_ = i_Pi - 1;

  if (p_ == materialPropertyTable_[indexPiMin_][0]) {
    indexPiPlus1Min_ = indexPiMin_;
    indexPiPlus1Max_ = indexPiMax_;
    f_p_ = 1;
  } else {
    indexPiPlus1Max_ = indexPiMin_ - 1;
    i_Pi = indexPiMin_ - 1;

    while
    (
        materialPropertyTable_[i_Pi][0] ==
        materialPropertyTable_[indexPiPlus1Max_][0]
    ) {
      i_Pi--;
      if (i_Pi == -1) {
        break;
      }
    }
    indexPiPlus1Min_ = i_Pi + 1;

    // f_p_=  (::log(p_) - ::log(materialPropertyTable_[indexPiMin_][0])) / (::log(materialPropertyTable_[indexPiPlus1Min_][0])- ::log(materialPropertyTable_[indexPiMin_][0]));
    f_p_ =
        (
            p_ -
            materialPropertyTable_[indexPiMin_][0]) /
        (
            materialPropertyTable_[indexPiPlus1Min_][0] -
            materialPropertyTable_[indexPiMin_][0]
        );
  }

  return f_p_;
}

// Returns M (and updates .. ) knowing p, T (interpolates in materialPropertyTable)
void Foam::LinearInterpolationMaterialPropertiesObject::update(scalar p_, scalar T_)
{
  T_temp_ = T_;
  p_temp_ = p_;

  // local variables
  scalar f_Ti_, cp_i_, cp_iPlus1_, h_i_, h_iPlus1_, ki_i_, ki_iPlus1_, kj_i_, kj_iPlus1_, kk_i_, kk_iPlus1_, eps_i_, eps_iPlus1_, alpha_i_, alpha_iPlus1_, hs_i_, hs_iPlus1_;

  // pressure indexes
  indexPi(p_);

  // temperature indexes for sub-table i
  f_Ti_ = indexTi(T_, indexPiMin_, indexPiMax_);

  // B'c and hw for sub-table i
  cp_i_ = f_Ti_ * materialPropertyTable_[indexTiPlus1_][2] + (1 - f_Ti_) * materialPropertyTable_[indexTi_][2];
  h_i_ = f_Ti_ * materialPropertyTable_[indexTiPlus1_][3] + (1 - f_Ti_) * materialPropertyTable_[indexTi_][3];
  ki_i_ = f_Ti_ * materialPropertyTable_[indexTiPlus1_][4] + (1 - f_Ti_) * materialPropertyTable_[indexTi_][4];
  kj_i_ = f_Ti_ * materialPropertyTable_[indexTiPlus1_][5] + (1 - f_Ti_) * materialPropertyTable_[indexTi_][5];
  kk_i_ = f_Ti_ * materialPropertyTable_[indexTiPlus1_][6] + (1 - f_Ti_) * materialPropertyTable_[indexTi_][6];
  eps_i_ = f_Ti_ * materialPropertyTable_[indexTiPlus1_][7] + (1 - f_Ti_) * materialPropertyTable_[indexTi_][7];
  alpha_i_ = f_Ti_ * materialPropertyTable_[indexTiPlus1_][8] + (1 - f_Ti_) * materialPropertyTable_[indexTi_][8];
  hs_i_ = f_Ti_ * materialPropertyTable_[indexTiPlus1_][9] + (1 - f_Ti_) * materialPropertyTable_[indexTi_][9];

  // temperature indexes for sub-table i+1
  f_Ti_ = indexTi(T_, indexPiPlus1Min_, indexPiPlus1Max_);

  // B'c for sub-table i,j+1
  cp_iPlus1_ =
      f_Ti_ * materialPropertyTable_[indexTiPlus1_][2] +
      (1 - f_Ti_) * materialPropertyTable_[indexTi_][2];
  h_iPlus1_ =
      f_Ti_ * materialPropertyTable_[indexTiPlus1_][3] +
      (1 - f_Ti_) * materialPropertyTable_[indexTi_][3];
  ki_iPlus1_ =
      f_Ti_ * materialPropertyTable_[indexTiPlus1_][4] +
      (1 - f_Ti_) * materialPropertyTable_[indexTi_][4];
  kj_iPlus1_ =
      f_Ti_ * materialPropertyTable_[indexTiPlus1_][5] +
      (1 - f_Ti_) * materialPropertyTable_[indexTi_][5];
  kk_iPlus1_ =
      f_Ti_ * materialPropertyTable_[indexTiPlus1_][6] +
      (1 - f_Ti_) * materialPropertyTable_[indexTi_][6];
  eps_iPlus1_ =
      f_Ti_ * materialPropertyTable_[indexTiPlus1_][7] +
      (1 - f_Ti_) * materialPropertyTable_[indexTi_][7];
  alpha_iPlus1_ =
      f_Ti_ * materialPropertyTable_[indexTiPlus1_][8] +
      (1 - f_Ti_) * materialPropertyTable_[indexTi_][8];
  hs_iPlus1_ =
      f_Ti_ * materialPropertyTable_[indexTiPlus1_][9] +
      (1 - f_Ti_) * materialPropertyTable_[indexTi_][9];

  // interpolations
  cp_    = (f_p_ * cp_iPlus1_ + (1 - f_p_) * cp_i_);
  h_     = (f_p_ * h_iPlus1_ + (1 - f_p_) * h_i_);
  ki_    = (f_p_ * ki_iPlus1_ + (1 - f_p_) * ki_i_);
  kj_    = (f_p_ * kj_iPlus1_ + (1 - f_p_) * kj_i_);
  kk_    = (f_p_ * kk_iPlus1_ + (1 - f_p_) * kk_i_);
  eps_   = (f_p_ * eps_iPlus1_ + (1 - f_p_) * eps_i_);
  alpha_ = (f_p_ * alpha_iPlus1_ + (1 - f_p_) * alpha_i_);
  hs_    = (f_p_ * hs_iPlus1_ + (1 - f_p_) * hs_i_);

}

Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::cp(scalar p_, scalar T_)
{
  if ( (T_temp_ != T_) || (p_temp_ != p_))
    update(p_, T_);

  return cp_;
}

Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::h(scalar p_, scalar T_)
{
  if ((T_temp_ != T_) || (p_temp_ != p_))
    update(p_, T_);

  return h_;
}

Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::ki(scalar p_, scalar T_)
{
  if ((T_temp_ != T_) || (p_temp_ != p_))
    update(p_, T_);

  return ki_;
}

Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::kj(scalar p_, scalar T_)
{
  if ((T_temp_ != T_) || (p_temp_ != p_))
    update(p_, T_);

  return kj_;
}

Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::kk(scalar p_, scalar T_)
{
  if ((T_temp_ != T_) || (p_temp_ != p_))
    update(p_, T_);

  return kk_;
}

Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::eps(scalar p_, scalar T_)
{
  if ((T_temp_ != T_) || (p_temp_ != p_))
    update(p_, T_);

  return eps_;
}

Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::alpha(scalar p_, scalar T_)
{
  if ((T_temp_ != T_) || (p_temp_ != p_))
    update(p_, T_);

  return alpha_;
}

Foam::scalar Foam::LinearInterpolationMaterialPropertiesObject::hs(scalar p_, scalar T_)
{
  if ((T_temp_ != T_) || (p_temp_ != p_))
    update(p_, T_);

  return hs_;
}
// ************************************************************************* //
