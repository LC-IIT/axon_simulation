TITLE hhy.mod   squid sodium, potassium, and leak channels
 
COMMENT
 This is the original Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application 
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set a squid-appropriate temperature
 (e.g. in HOC: "celsius=6.3" or in Python: "h.celsius=6.3").
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
 
 Slow inactivation is modeled by s-gate as in Zang et al. Curr Biol 2023:
 DOI: 10.1016/j.cub.2023.03.043, https://modeldb.science/267610
 and one additional parameter a0s_fold, set to 0.3
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
:		(S) = (siemens)
}

NEURON {
        SUFFIX hhy

        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk,a2,s_init,smin,tslow,a0s_fold

        GLOBAL minf, hinf, ninf, mtau, htau, ntau
:	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .12 (S/cm2)	<0,1e9>
        gkbar = .036 (S/cm2)	<0,1e9>
        gl = .0003 (S/cm2)	<0,1e9>
        el = -54.3 (mV)
        
        gamma = 0.0000025 (/ms)
        delta = 0.25 (/ms)
        
    	vhalfs=-60	(mV)		: slow inact.
        a0s=0.0003	(ms)		: a0s=b0s
        smin=10		(ms)
        vvh=-58		(mV) 
        vvs=2		(mV)
        a2=1		(1)		: 1=no inact., 0=max inact.
        s_init = 1
        tslow = 0
        a0s_fold = 1                    : scaling factor for the slow recovery tau 
}
 
STATE {
        m h n s
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

		gna (S/cm2)
		gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf sinf
		mtau (ms)
		htau (ms)
		ntau (ms)
		stau (ms)
}

: currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h*s
		ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
		ik = gk*(v - ek)      
        il = gl*(v - el)
}

: states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
		s' = (sinf-s)/stau

}

UNITSOFF

INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	s = s_init       
}

PROCEDURE rates(v(mV)) { 

        LOCAL  alpha, beta, sum, q10,cs,avail

        q10 = 3^((celsius - 6.3)/10)
                :"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10)
        beta =  4 * exp(-(v+65)/18)
        sum = alpha + beta
		mtau = 1/(q10*sum)
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20)
        beta = 1 / (exp(-(v+35)/10) + 1)
        sum = alpha + beta
		htau = 1/(q10*sum)
        hinf = alpha/sum
                :"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
		sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
 		
        cs = 1/(1+exp((v-vvh)/vvs))
        avail = 1
        if (t>2000) {avail = a2}
        sinf = cs+avail*(1-cs)
        if (t<2000) {sinf = s_init}
        
        alpha = exp(0.45*(v-vhalfs))
        beta = exp(0.09*(v-vhalfs))
        stau = beta/(a0s_fold*a0s*(1+alpha))
        
        if (stau<smin) {stau=smin}
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON