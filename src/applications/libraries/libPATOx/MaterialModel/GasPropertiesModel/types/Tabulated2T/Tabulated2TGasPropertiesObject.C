#include "Tabulated2TGasPropertiesObject.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Tabulated2TGasPropertiesObject::Tabulated2TGasPropertiesObject
(
    RectangularMatrix<scalar>& gasPropertyTable
)
  :
gasPropertyTable_(gasPropertyTable),
indexPiMin_(0),
indexPiMax_(0),
indexPiPlus1Min_(0),
indexPiPlus1Max_(0),
f_p_(0.0),
indexTi_(0.0),
indexTiPlus1_(0.0),
f_T_(0.0),
M_(0.0),
cp_(0.0),
//gamma_(0.0),
h_(0.0),
mu_(0.0),
//rho_(0.0),
k_(0.0),
T_temp_(0.0),
p_temp_(0.0),
small_(0.0)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Tabulated2TGasPropertiesObject::~Tabulated2TGasPropertiesObject()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// indexTi: Returns the Temperature interpolation factor and updates indexTi_ and indexTiPlus1_
Foam::scalar Foam::Tabulated2TGasPropertiesObject::indexTi(scalar T_, int i_min_, int i_max_)
{

  if
  (
      (T_ - gasPropertyTable_[i_min_][1]) < -small_
      || (T_ - gasPropertyTable_[i_max_][1]) > small_
  ) {
//    Info << "Temperature out of gasPropertyTable range" << nl;

    if ( (T_ - gasPropertyTable_[i_min_][1]) < - small_ ) {
//      Info << "T = " << T_
//           << " < " << gasPropertyTable_[i_min_][1]
//           << nl;
      T_ = gasPropertyTable_[i_min_][1];
    }

    if ( (T_ - gasPropertyTable_[i_max_][1]) > small_ ) {
//      Info << "T = " << T_
//           << " > " << gasPropertyTable_[i_max_][1]
//           << nl;
      T_ =      gasPropertyTable_[i_max_][1];
    }
  }

  int i_Ti(i_min_);

  while (gasPropertyTable_[i_Ti][1] < T_) {
    i_Ti++;
  }

  indexTi_ = i_Ti - 1;
  if (indexTi_<0) {
    indexTi_=0;
  }
  if (T_ == gasPropertyTable_[indexTi_][1]) {
    indexTiPlus1_ = indexTi_;
    f_T_ = 1;
    return f_T_;
  } else {
    indexTiPlus1_ = indexTi_ + 1;
    f_T_ =
        (
            T_ -
            gasPropertyTable_[indexTi_][1]
        ) /
        (
            gasPropertyTable_[indexTiPlus1_][1] -
            gasPropertyTable_[indexTi_][1]
        );
  }

  // Info << "f_T = " << f_T_ << nl;
  return f_T_;
}

// indexPi: Returns the interpolation factor f and updates IndexPiMin_, IndexPiMax_, IndexPiPlus1Min_, IndexPiPlus1Max_
Foam::scalar Foam::Tabulated2TGasPropertiesObject::indexPi(scalar p_)
{
  int i_Pi(0);

  if
  (
      (p_ - gasPropertyTable_[0][0]) > small_
      || (p_ - gasPropertyTable_[gasPropertyTable_.m() - 1][0]) < -small_
  ) {
//    Info << "Pressure out of gasPropertyTable range" << nl;

    if (p_ - gasPropertyTable_[0][0] > small_) {
//      Info << "p = " << p_
//           << " > " << gasPropertyTable_[0][0]
//           << nl;
      p_ = gasPropertyTable_[0][0];
    }

    if (p_ - gasPropertyTable_[gasPropertyTable_.m() - 1][0] < -small_) {
//      Info << "p = " << p_
//           << " < " << gasPropertyTable_[gasPropertyTable_.m() - 1][0]
//           << nl;
      p_ = gasPropertyTable_[gasPropertyTable_.m() - 1][0];
    }

  }

  while (gasPropertyTable_[i_Pi][0] > p_) {
    i_Pi++;
  }

  indexPiMin_ = i_Pi;

  while (gasPropertyTable_[i_Pi][0] == gasPropertyTable_[indexPiMin_][0]) {
    i_Pi++;
    if (i_Pi == gasPropertyTable_.m()) {
      break;
    }
  }

  indexPiMax_ = i_Pi - 1;
  if (indexPiMax_<0) {
    indexPiMax_=0;
  }
  if (p_ == gasPropertyTable_[indexPiMin_][0]) {
    indexPiPlus1Min_ = indexPiMin_;
    indexPiPlus1Max_ = indexPiMax_;
    f_p_ = 1;
  } else {
    indexPiPlus1Max_ = indexPiMin_ - 1;
    if (indexPiPlus1Max_<0) {
      indexPiPlus1Max_=0;
    }
    i_Pi = indexPiMin_ - 1;

    while
    (
        gasPropertyTable_[i_Pi][0]
        == gasPropertyTable_[indexPiPlus1Max_][0]
    ) {
      i_Pi--;

      if (i_Pi == -1)
        break;
    }
    indexPiPlus1Min_ = i_Pi + 1;
    if (indexPiPlus1Min_<0) {
      indexPiPlus1Min_=0;
    }
//    f_p_=  (::log(p_) - ::log(gasPropertyTable_[indexPiMin_][0])) / (::log(gasPropertyTable_[indexPiPlus1Min_][0])- ::log(gasPropertyTable_[indexPiMin_][0]));

    f_p_ =
        (
            p_ -
            gasPropertyTable_[indexPiMin_][0]
        ) /
        (
            gasPropertyTable_[indexPiPlus1Min_][0] -
            gasPropertyTable_[indexPiMin_][0]
        );
  }

  return f_p_;
}

// Returns M (and updates .. ) knowing p, T (interpolates in gasPropertyTable)
Foam::scalar Foam::Tabulated2TGasPropertiesObject::M(scalar p_, scalar T_)
{
  T_temp_ = T_;
  p_temp_ = p_;
  // local variables
  scalar f_Ti_, M_i_, M_iPlus1_, cp_i_, cp_iPlus1_, h_i_, h_iPlus1_, mu_i_, mu_iPlus1_, k_i_, k_iPlus1_;
  //scalar cp_i_, cp_iPlus1_, gamma_i_, gamma_iPlus1_, rho_i_, rho_iPlus1_

  // pressure indexes
  indexPi(p_);

  // temperature indexes for sub-table i
  f_Ti_ = indexTi(T_, indexPiMin_, indexPiMax_);

  // B'c and hw for sub-table i
  M_i_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][2] + (1 - f_Ti_) * gasPropertyTable_[indexTi_][2];
  cp_i_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][3] + (1-f_Ti_) * gasPropertyTable_[indexTi_][3];
  //gamma_i_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][4] + (1-f_Ti_) * gasPropertyTable_[indexTi_][4];
  h_i_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][4] + (1 - f_Ti_) * gasPropertyTable_[indexTi_][4];
  mu_i_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][5] + (1 - f_Ti_) * gasPropertyTable_[indexTi_][5];
  //rho_i_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][7] + (1-f_Ti_) * gasPropertyTable_[indexTi_][6];
  k_i_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][6] + (1-f_Ti_) * gasPropertyTable_[indexTi_][6];

  // temperature indexes for sub-table i+1
  f_Ti_ = indexTi(T_, indexPiPlus1Min_, indexPiPlus1Max_);

  // B'c for sub-table i,j+1
  M_iPlus1_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][2] + (1 - f_Ti_) * gasPropertyTable_[indexTi_][2];
  cp_iPlus1_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][3] + (1-f_Ti_) * gasPropertyTable_[indexTi_][3];
  //gamma_iPlus1_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][4] + (1-f_Ti_) * gasPropertyTable_[indexTi_][4];
  h_iPlus1_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][4] + (1 - f_Ti_) * gasPropertyTable_[indexTi_][4];
  mu_iPlus1_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][5] + (1 - f_Ti_) * gasPropertyTable_[indexTi_][5];
  //rho_iPlus1_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][7] + (1-f_Ti_) * gasPropertyTable_[indexTi_][6];
  k_iPlus1_ = f_Ti_ * gasPropertyTable_[indexTiPlus1_][6] + (1-f_Ti_) * gasPropertyTable_[indexTi_][6];

  M_ =  (f_p_ * M_iPlus1_ + (1 - f_p_) * M_i_);
  cp_=  (f_p_ * cp_iPlus1_ + (1-f_p_) * cp_i_);
  //gamma_=  (f_p_ * gamma_iPlus1_ + (1-f_p_) * gamma_i_);
  h_ =  (f_p_ * h_iPlus1_ + (1 - f_p_) * h_i_);
  mu_ =  (f_p_ * mu_iPlus1_ + (1 - f_p_) * mu_i_);
  //rho_=  (f_p_ * rho_iPlus1_ + (1-f_p_) * rho_i_);
  k_=  (f_p_ * k_iPlus1_ + (1-f_p_) * k_i_);

  return M_;
}

Foam::scalar Foam::Tabulated2TGasPropertiesObject::cp(scalar p_, scalar T_)
{
  if
  (
      (T_temp_!=T_)
      || (p_temp_!=p_)
  ) {
    M(p_,T_);
  }

  return cp_;
}

Foam::scalar Foam::Tabulated2TGasPropertiesObject::h(scalar p_, scalar T_)
{
  if
  (
      (T_temp_ != T_)
      || (p_temp_ != p_)
  ) {
    M(p_, T_);
  }

  return h_;
}

Foam::scalar Foam::Tabulated2TGasPropertiesObject::mu(scalar p_, scalar T_)
{
  if
  (
      (T_temp_ != T_)
      || (p_temp_ != p_)
  ) {
    M(p_, T_);
  }

  return mu_;
}

Foam::scalar Foam::Tabulated2TGasPropertiesObject::k(scalar p_, scalar T_)
{
  if
  (
      (T_temp_!=T_)
      || (p_temp_!=p_)
  ) {
    M(p_,T_);
  }

  return k_;
}

/*
Foam::scalar Foam::Tabulated2TGasPropertiesObject::gamma(scalar p_, scalar T_)
{
    if ((T_temp_!=T_) || (p_temp_!=p_)){M(p_,T_);}

    return gamma_;
}

Foam::scalar Foam::Tabulated2TGasPropertiesObject::rho(scalar p_, scalar T_)
{
    if
    (
	(T_temp_!=T_)
	|| (p_temp_!=p_)
    ) {
      M(p_,T_);
    }

    return rho_;
}
*/

// ************************************************************************* //
