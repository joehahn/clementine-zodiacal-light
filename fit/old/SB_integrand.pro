function SB_integrand, scat_angle, PSI=psi, H=h, BETA2=beta2

;Calculate the integrand of the zodiacal-light
;surface-brightness integral
;psi(scat_angle)*h(beta)*sin(scat_angle)^nu,
;where psi(scat_angle)=phase-law, h(beta)=latitude distribution,
;and scat_angle=scattering angle.

common share, beta, sigma
common ZL, phi, theta, epsilon, nu, err

;Hong's phase-function psi
psi=phase_law(scat_angle,nu)

;solar latitude beta'
beta=epsilon
if (epsilon gt 1.0e-6) then $
  beta=abs(asin(sin(scat_angle-epsilon)*sin(theta)/sin(epsilon)))
if (beta le 0.0) then beta=0.0
beta2=beta

;dust latitude distribution h(beta) normalized
;such that h(beta=0)=1
Po2=1.57080
h =qromo('latitude_integrand',beta,Po2,double=0,eps=err)
beta=0.0
h0=qromo('latitude_integrand',beta,Po2,double=0,eps=err)
h=h/h0

;surface-brightness integrand
SB=psi*h*(sin(scat_angle)^nu)

return,SB
end

