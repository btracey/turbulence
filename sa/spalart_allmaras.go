package sa

import "math"

/*
type SAConstantList str
	CB1                 float64
	Sigma               float64
	CB2                 float64
	Kappa               float64
	CW2                 float64
	CW3                 float64
	CV1                 float64
	CT3                 float64
	CT4                 float64
	CW1                 float64
	CV1Cubed            float64
	CW3Sixthed          float64
	CB1OverKappaSquared float64
	CB2OverSigma        float64
}

*/

const (
	CB1                 = 0.1355
	Sigma               = 2.0 / 3.0
	CB2                 = 0.622
	Kappa               = 0.41
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

type SA struct {
	Unlimited bool // Should caps be put on the variables
	// Inputs
	NDim                int
	Nu                  float64
	NuHat               float64
	DNuHatDX            []float64
	DUIdXJ              [][]float64
	WallDistance        float64
	LimitedWallDistance float64

	// Intermediate quantities
	W     [][]float64
	Chi   float64
	Ft2   float64
	Omega float64
	SHat  float64
	Fv1   float64
	Fv2   float64

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

func (s *SA) computeFv1() {
	chicubed := math.Pow(s.Chi, 3)
	s.Fv1 = chicubed / (chicubed + CV1Cubed)
}

func (s *SA) computeFt2() {
	s.Ft2 = CT3 * math.Exp(-CT4*s.Chi*s.Chi)
}

func (s *SA) computeFv2() {
	s.computeFv1()
	s.Fv2 = 1 - s.Chi/(1+s.Chi*s.Fv1)
}

func (s *SA) computeSHat() {
	s.computeOmega()
	s.computeFv2()
	s.SHat = s.Omega + (s.NuHat/(Kappa*Kappa*s.WallDistance*s.WallDistance))*s.Fv2
	//if s.SHat < 1e-8 {
	//	s.SHat = 1e-8
	//}
}

func (s *SA) computeOmega() {
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

func (s *SA) computeR() {
	r := s.NuHat / (s.SHat * Kappa * Kappa * s.WallDistance * s.WallDistance)
	if r > 10 {
		r = 10
	}
	s.R = r
}

func (s *SA) computeG() {
	s.ComputeR()
	s.G = s.R + CW2*(math.Pow(s.R, 6)-s.R)
}

func (s *SA) computeFw() {
	s.ComputeG()
	limiter := (1 + CW3Sixthed) / (math.Pow(s.G, 6) + CW3Sixthed)
	s.Fw = s.G * math.Pow(limiter, 1.0/6.0)
}

func (s *SA) computeProduction() {
	s.Production = CB1 * (1 - s.Ft2) * s.SHat * s.NuHat
}

func (s *SA) computeDestruction() {
	s.Production = (CW1*s.Fw - Cb1OverKappaSquared*s.Ft2) * (s.NuHat / s.WallDistance) * (s.NuHat / s.WallDistance)
}

func (s *SA) computeCrossProduction() {
	s.CrossProduction = 0
	for i := 0; i < s.NDim; i++ {
		s.CrossProduction += Cb2OverSigma * s.DNuHatDX[i] * s.DNuHatDX[i]
	}
}

// Compute all of the values
func (s *SA) Source() float64 {
	s.LimitedWallDistance = s.WallDistance
	if !s.Unlimited {
		if s.LimitedWallDistance < 1e-10 {
			s.LimitedWallDistance = 1e-10
		}
	}
	s.Chi = s.NuHat / s.Nu
	s.computeFt2()
	s.computeSHat()
	if !s.Unlimited {
		if s.SHat < 1e-10 {
			s.SHat = 1e-10
		}
	}
	s.computeFw()
	s.computeFv1()
	s.computeProduction()
	s.computeDestruction()
	s.computeCrossProduction()
	s.SourceTerm = s.Production - s.Destruction + s.CrossDiffusion
	return s.SourceTerm
}
