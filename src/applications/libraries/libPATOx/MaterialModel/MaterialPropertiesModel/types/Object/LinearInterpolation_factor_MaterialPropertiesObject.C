#include "LinearInterpolation_factor_MaterialPropertiesObject.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LinearInterpolation_Factor_MaterialPropertiesObject::LinearInterpolation_Factor_MaterialPropertiesObject
(
    const LinearInterpolationMaterialPropertiesObject& obj,
    List<scalar>& factors
)
  :
LinearInterpolationMaterialPropertiesObject(obj),
factors_(factors)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LinearInterpolation_Factor_MaterialPropertiesObject::~LinearInterpolation_Factor_MaterialPropertiesObject()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Returns M (and updates .. ) knowing p, T (interpolates in materialPropertyTable)
void Foam::LinearInterpolation_Factor_MaterialPropertiesObject::update(scalar p_, scalar T_)
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
  cp_    = (f_p_ * cp_iPlus1_ + (1 - f_p_) * cp_i_)*factors_[0];
  h_     = (f_p_ * h_iPlus1_ + (1 - f_p_) * h_i_)*factors_[1];
  ki_    = (f_p_ * ki_iPlus1_ + (1 - f_p_) * ki_i_)*factors_[2];
  kj_    = (f_p_ * kj_iPlus1_ + (1 - f_p_) * kj_i_)*factors_[3];
  kk_    = (f_p_ * kk_iPlus1_ + (1 - f_p_) * kk_i_)*factors_[4];
  eps_   = (f_p_ * eps_iPlus1_ + (1 - f_p_) * eps_i_)*factors_[5];
  alpha_ = (f_p_ * alpha_iPlus1_ + (1 - f_p_) * alpha_i_)*factors_[6];
  hs_    = (f_p_ * hs_iPlus1_ + (1 - f_p_) * hs_i_);

}

// ************************************************************************* //
