  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    /* Standard CDM: rho_cdm = Omega0_cdm * H0^2 / a^3 */
    double rho_cdm_std = pba->Omega0_cdm * pow(pba->H0,2) / pow(a,3);
    
    /* Holographic background modification - conservative approach */
    if (pba->interaction_beta != 0 && a > 0.1) {
      /* Small late-time enhancement from DE->CDM energy transfer */
      double beta = pba->interaction_beta;
      /* Smooth turn-on after a=0.1, max effect at a=1 */
      double x = (a - 0.1) / 0.9;  /* 0 to 1 for a in [0.1, 1] */
      if (x > 1.0) x = 1.0;
      if (x < 0.0) x = 0.0;
      /* Small enhancement: 1 + 0.05*beta at a=1 */
      double enhancement = 1.0 + 0.05 * beta * x * x;
      pvecback[pba->index_bg_rho_cdm] = rho_cdm_std * enhancement;
    } else {
      pvecback[pba->index_bg_rho_cdm] = rho_cdm_std;
    }
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_cdm];
  }
