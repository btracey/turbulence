package sa

import "math"

const (
	CB1                 = 0.1355
	Sigma               = 2.0 / 3.0
	CB2                 = 0.622
	Kappa               = 0.41
	K2                  = Kappa * Kappa
	CW2                 = 0.3
	CW3                 = 2
	CV1                 = 7.1
	CV1Cubed            = CV1 * CV1 * CV1
	CW3Sixthed          = CW3 * CW3 * CW3 * CW3 * CW3 * CW3
	Cb1OverKappaSquared = CB1 / (Kappa * Kappa)
	CT3                 = 1.2
	CT4                 = 0.5
	CW1                 = 0.1355/(0.41*0.41) + (1+0.622)/(2.0/3.0)
	CB2OverSigma        = CB2 / Sigma
)

// TODO: At some point, it should be SA.Compute(). Also, should replace with a
// 2-D and 3-D

type SA struct {
	Unlimited bool // Should caps be put on the variables
	// Inputs
	NDim         int
	Nu           float64
	NuHat        float64
	DNuHatDX     []float64
	DUIdXJ       [][]float64
	WallDistance float64

	// Intermediate quantities
	W                   [][]float64
	Chi                 float64
	Ft2                 float64
	Omega               float64
	SHat                float64
	InvSHat             float64
	Fv1                 float64
	Fv2                 float64
	LimitedWallDistance float64

	Production      float64
	Destruction     float64
	CrossProduction float64
	SourceTerm      float64

	G  float64
	Fw float64
	R  float64

	ThirdTerm     float64
	MatDerivDNuDt float64
}

func (s *SA) ComputeFv1() {
	chicubed := math.Pow(s.Chi, 3)
	s.Fv1 = chicubed / (chicubed + CV1Cubed)
}

func (s *SA) ComputeFt2() {
	s.Ft2 = CT3 * math.Exp(-CT4*s.Chi*s.Chi)
}

func (s *SA) ComputeFv2() {
	s.ComputeFv1()
	s.Fv2 = 1 - s.Chi/(1+s.Chi*s.Fv1)
}

func (s *SA) ComputeSHat() {
	s.ComputeOmega()
	s.ComputeFv2()
	d2 := s.WallDistance * s.LimitedWallDistance
	invk2d2 := 1.0 / (K2 * d2)
	s.SHat = s.Omega + s.NuHat*s.Fv2*invk2d2
}

func (s *SA) ComputeOmega() {
	// First, compute W
	for i := range s.W {
		for j := range s.W[i] {
			s.W[i][j] = 0.5 * (s.DUIdXJ[i][j] - s.DUIdXJ[j][i])
		}
	}

	// Omega = sqrt( 2 Wij Wij)
	for i := range s.W {
		for j := range s.W[i] {
			s.Omega += s.W[i][j] * s.W[i][j]
		}
	}
	s.Omega *= 2
	s.Omega = math.Sqrt(s.Omega)
}

func (s *SA) ComputeR() {
	r := s.NuHat * s.InvSHat / (Kappa * Kappa * s.LimitedWallDistance * s.LimitedWallDistance)
	if r > 10 {
		r = 10
	}
	s.R = r
}

func (s *SA) ComputeG() {
	s.ComputeR()
	s.G = s.R + CW2*(math.Pow(s.R, 6)-s.R)
}

func (s *SA) ComputeFw() {
	s.ComputeG()
	limiter := (1 + CW3Sixthed) / (math.Pow(s.G, 6) + CW3Sixthed)
	s.Fw = s.G * math.Pow(limiter, 1.0/6.0)
}

func (s *SA) ComputeProduction() {
	//s.Production = CB1 * (1 - s.Ft2) * s.SHat * s.NuHat
	s.Production = CB1 * s.SHat * s.NuHat
}

func (s *SA) ComputeDestruction() {
	//s.Destruction = (CW1*s.Fw - Cb1OverKappaSquared*s.Ft2) * (s.NuHat / s.WallDistance) * (s.NuHat / s.WallDistance)
	s.Destruction = CW1 * s.Fw * (s.NuHat / s.LimitedWallDistance) * (s.NuHat / s.LimitedWallDistance)
}

func (s *SA) ComputeCrossProduction() {
	s.CrossProduction = 0
	for i := 0; i < s.NDim; i++ {
		s.CrossProduction += CB2OverSigma * s.DNuHatDX[i] * s.DNuHatDX[i]
	}
}

// Compute all of the values
func (s *SA) Source() float64 {
	if s.WallDistance == 0 {
		return 0
	}
	s.LimitedWallDistance = s.WallDistance
	if !s.Unlimited {
		if s.LimitedWallDistance < 1e-10 {
			s.LimitedWallDistance = 1e-10
		}
	}

	if len(s.DNuHatDX) != s.NDim {
		panic("nu deriv length mismatch")
	}
	if len(s.DUIdXJ) != s.NDim {
		panic("velocity deriv length mismatch")
	}
	for _, v := range s.DUIdXJ {
		if len(v) != s.NDim {
			panic("velocity deriv length mismatch")
		}
	}

	s.W = make([][]float64, s.NDim)
	for i := range s.W {
		s.W[i] = make([]float64, s.NDim)
	}
	s.Chi = s.NuHat / s.Nu
	s.ComputeFt2()
	s.ComputeSHat()

	s.InvSHat = 1.0 / math.Max(s.SHat, 10e-10)
	s.ComputeFw()
	s.ComputeFv1()
	s.ComputeProduction()
	s.ComputeDestruction()
	s.ComputeCrossProduction()
	s.SourceTerm = s.Production - s.Destruction + s.CrossProduction
	return s.SourceTerm
}

const (
	//NuAir = 1.48e-5
	NuAir = 1.0
)

/*
func NuAlt(nu float64) float64 {
	return NuAir
}
*/
func NuAltScale(nu float64) float64 {
	return nu / NuAir
}

func NuHatAlt(nu, nuhat float64) float64 {
	nuAlt := NuAltScale(nu)
	return nuhat / nuAlt
}

func NuHatAltStar(nu, nuhat float64) float64 {
	nuHatAlt := NuHatAlt(nu, nuhat)
	return nuHatAlt + 3*NuAir
}

func DistAlt(dist, nu float64) float64 {
	//nuAlt := NuAltScale(nu)
	//return dist / nuAlt
	return dist
}

func DistAltStar(dist, nu float64) float64 {
	//distAlt := DistAlt(dist, nu)
	//return math.Min(distAlt, 300*NuAir) // 100 * nuhat_0

	return math.Min(dist, 300*nu)
}

func SourceAlt(source, nu float64) float64 {
	nutilde := NuAltScale(nu)
	return source / (nutilde * nutilde)
}

func SourceNondimerAlt(nu, nuhat, dist float64) float64 {
	nuhatAlt := NuHatAlt(nu, nuhat)
	distAlt := DistAlt(dist, nu)
	return (nuhatAlt / distAlt) * (nuhatAlt / distAlt)
}

func SourceNondimerAltStar(nu, nuhat, dist float64) float64 {
	nuhatAltStar := NuHatAltStar(nu, nuhat)
	distStar := DistAltStar(dist, nu)
	sourceNondimer := (nuhatAltStar / distStar) * (nuhatAltStar / distStar)
	/*
		if sourceNondimer > 500000 && dist != 0 {
			fmt.Println("nu", nu)
			fmt.Println("nuhat", nuhat)
			fmt.Println("dist", dist)
			fmt.Println("nuhatalt", NuHatAlt(nu, nuhat))
			fmt.Println("nu hat alt Star", nuhatAltStar)
			fmt.Println(distStar)
			fmt.Println(sourceNondimer)
			os.Exit(1)
		}
	*/
	return sourceNondimer
}

func InvSourceNondimerAltStar(nu, nuhat, dist float64) float64 {
	return 1 / SourceNondimerAltStar(nu, nuhat, dist)
}

func NondimSourceAltStar(source, nu, nuhat, dist float64) float64 {
	sourcenondimer := SourceNondimerAltStar(nu, nuhat, dist)
	sourceAlt := SourceAlt(source, nu)
	return sourceAlt / sourcenondimer
}

func NondimSourceAlt(source, nu, nuhat, dist float64) float64 {
	sourcenondimer := SourceNondimerAlt(nu, nuhat, dist)
	sourceAlt := SourceAlt(source, nu)
	return sourceAlt / sourcenondimer
}

func OmegaAlt(omega, nu float64) float64 {
	nutilde := NuAltScale(nu)
	return omega / nutilde
}

func OmegaAltNondimRatio(nu, nuhat, dist float64) float64 {
	nuhatAltStar := NuHatAltStar(nu, nuhat)
	distStar := DistAltStar(dist, nu)
	ratOrig := nuhat / (dist * dist)
	ratNew := nuhatAltStar / (distStar * distStar)
	return ratNew / ratOrig
}

func LogOmegaAltNondimRatio(nu, nuhat, dist float64) float64 {
	nuhatAltStar := NuHatAltStar(nu, nuhat)
	distStar := DistAltStar(dist, nu)
	ratOrig := nuhat / (dist * dist)
	ratNew := nuhatAltStar / (distStar * distStar)
	return math.Log(ratNew) - math.Log(ratOrig)
}

func SourceAltNondimRatio(nu, nuhat, dist float64) float64 {
	nuhatAltStar := NuHatAltStar(nu, nuhat)
	distStar := DistAltStar(dist, nu)
	ratOrig := (nuhat / dist) * (nuhat / dist)
	ratNew := (nuhatAltStar / distStar) * (nuhatAltStar / distStar)
	return ratNew / ratOrig
}

func LogSourceAltNondimRatio(nu, nuhat, dist float64) float64 {
	nuhatAltStar := NuHatAltStar(nu, nuhat)
	distStar := DistAltStar(dist, nu)
	ratOrig := (nuhat / dist) * (nuhat / dist)
	ratNew := (nuhatAltStar / distStar) * (nuhatAltStar / distStar)
	return math.Log(ratNew) - math.Log(ratOrig)
}

func OmegaNondimerAlt(nu, nuhat, dist float64) float64 {
	nuhatAlt := NuHatAlt(nu, nuhat)
	distAlt := DistAlt(dist, nu)
	return (1 / distAlt) * (nuhatAlt / distAlt)
}

func OmegaNondimerAltStar(nu, nuhat, dist float64) float64 {
	nuhatAltStar := NuHatAltStar(nu, nuhat)
	distStar := DistAltStar(dist, nu)
	return (1 / distStar) * (nuhatAltStar / distStar)
}

func InvOmegaNondimerAltStar(nu, nuhat, dist float64) float64 {
	return 1 / OmegaNondimerAltStar(nu, nuhat, dist)
}

func NondimOmegaAlt(omega, nu, nuhat, dist float64) float64 {
	omegaNondimer := OmegaNondimerAlt(nu, nuhat, dist)
	omegaAlt := OmegaAlt(omega, nu)
	return omegaAlt / omegaNondimer
}

func NondimOmegaAltStar(omega, nu, nuhat, dist float64) float64 {
	omegaNondimer := OmegaNondimerAltStar(nu, nuhat, dist)
	omegaAlt := OmegaAlt(omega, nu)
	return omegaAlt / omegaNondimer
}

func ChiAlt(chi float64) float64 {
	return math.Min(chi, 100)
}

func NuGradMagAlt(nuGradMagAlt, nu float64) float64 {
	nutilde := NuAltScale(nu)
	return nuGradMagAlt / (nutilde * nutilde)
}

func NondimNuGradMagAlt(grad, nu, nuhat, dist float64) float64 {
	sourcenondimer := SourceNondimerAlt(nu, nuhat, dist)
	nuGradAlt := NuGradMagAlt(grad, nu)
	return nuGradAlt / sourcenondimer
}

func NondimNuGradMagAltStar(grad, nu, nuhat, dist float64) float64 {
	sourcenondimer := SourceNondimerAltStar(nu, nuhat, dist)
	nuGradAlt := NuGradMagAlt(grad, nu)

	/*
		fmt.Println(nuGradAlt)

		if nuGradAlt/sourcenondimer > 170 {
			fmt.Println()
			fmt.Println("dist", dist)
			fmt.Println("dist alt", DistAlt(dist, nu))
			fmt.Println("dist alt *", DistAltStar(dist, nu))
			fmt.Println("nu", nu)
			fmt.Println("nuscale", NuAltScale(nu))
			fmt.Println("nuhat", nuhat)
			fmt.Println("nuhat alt *", NuHatAlt(nu, nuhat))
			fmt.Println("nuhat alt star ", NuHatAltStar(nu, nuhat))
			fmt.Println(sourcenondimer)
			fmt.Println(grad)
			fmt.Println(nuGradAlt)
			fmt.Println(nuGradAlt / sourcenondimer)
			os.Exit(1)
		}
	*/
	return nuGradAlt / sourcenondimer
}

func Fv1(chi float64) float64 {
	return chi * chi * chi / (chi*chi*chi + CV1Cubed)
}

func Fv2(chi, fv1 float64) float64 {
	return 1 - chi/(1+chi*fv1)
}

func Fw(g float64) float64 {
	limiter := (1 + CW3Sixthed) / (math.Pow(g, 6) + CW3Sixthed)
	return g * math.Pow(limiter, 1.0/6.0)
}

func G(r float64) float64 {
	return r + CW2*(math.Pow(r, 6)-r)
}
