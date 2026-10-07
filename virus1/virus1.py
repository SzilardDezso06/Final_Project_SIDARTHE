'''
S - Susceptible
I - Infected
D - Diagnosed
A - Ailing
R - Recognized
T - Threatened
H - Heded
E - Extinct

The parameters alpha, beta, gamma, and delta represent the disease transmission rates between a susceptible person and individuals who are infected, diagnosed, ailing, or recognized, respectively
(Generally, transmission is highest from infected individuals (alpha), lower from those showing symptoms (gamma), and lowest from those who are diagnosed or isolated (beta and delta). These rates 
can be actively reduced through social-distancing policies like lockdowns and remote work, while the risk from patients in proper intensive care units is considered negligible.)


The parameters epsilon and theta represent the detection rates for asymptomatic and symptomatic cases, respectively. 
(These values reflect the overall testing volume and contact-tracing efforts, with theta typically being higher because symptomatic individuals 
are more likely to seek testing. Public health officials can increase both rates by launching massive testing and tracing campaigns.)

The parameters zeta and eta represent the rates at which unaware and aware infected individuals, respectively, develop noticeable clinical symptoms. 
(While these values are normally comparable and depend heavily on the specific disease, they can be reduced over time through improved medical therapies and the population's acquisition of immunity.)

The parameters mu and nu represent the rates at which undetected and detected infected subjects, respectively, develop life-threatening symptoms. 
(Their values are highly dependent on available treatments (which can lower nu) and diagnosis patterns (since more severe cases are more likely to be detected, raising nu). 
Similar to the symptom-development rates, these parameters can be reduced through improved medical therapies and immunity.)

The parameter tau represents the mortality rate for infected individuals who have developed life-threatening symptoms. 
(This rate depends heavily on medical capabilities and can be reduced over time through improved therapies.)

The parameters llambda, kappa, xi, rho, and sigma represent the recovery rates for the five different classes of infected individuals. 
(These rates are generally comparable unless specific, effective treatments are introduced for diagnosed patients, which would cause certain rates to differ significantly. 
Like previous metrics, they can be increased through improved medical treatments and population immunity.)

'''

dS/dt = -S*(alpha*I + beta*D + gamma*A + delta*R)

dI/dt = S*(alpha*I + beta*D + gamma*A + delta*R) - (epsilon + zeta + llambda)*I

dD/dt = epsilon*I - (eta + rho)*D

dA/dt = zeta*I - (theta + mu + kappa)*A

dR/dt = eta*D + theta*A - (nu + xi)*R

dT/dt = mu*A + nu*R - (sigma + tau)*T

dH/dt = llambda*I + rho*D + kappa*A + xi*R + sigma*T

dE/dt = tau*T
