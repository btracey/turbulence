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
	d2 := s.WallDistance * s.LimitedWallDistance
	invk2d2 := 1.0 / (K2 * d2)
	s.SHat = s.Omega + s.NuHat*s.Fv2*invk2d2
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
	r := s.NuHat * s.InvSHat / (Kappa * Kappa * s.LimitedWallDistance * s.LimitedWallDistance)
	if r > 10 {
		r = 10
	}
	s.R = r
}

func (s *SA) computeG() {
	s.computeR()
	s.G = s.R + CW2*(math.Pow(s.R, 6)-s.R)
}

func (s *SA) computeFw() {
	s.computeG()
	limiter := (1 + CW3Sixthed) / (math.Pow(s.G, 6) + CW3Sixthed)
	s.Fw = s.G * math.Pow(limiter, 1.0/6.0)
}

func (s *SA) computeProduction() {
	//s.Production = CB1 * (1 - s.Ft2) * s.SHat * s.NuHat
	s.Production = CB1 * s.SHat * s.NuHat
}

func (s *SA) computeDestruction() {
	//s.Destruction = (CW1*s.Fw - Cb1OverKappaSquared*s.Ft2) * (s.NuHat / s.WallDistance) * (s.NuHat / s.WallDistance)
	s.Destruction = CW1 * s.Fw * (s.NuHat / s.LimitedWallDistance) * (s.NuHat / s.LimitedWallDistance)
}

func (s *SA) computeCrossProduction() {
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
	s.computeFt2()
	s.computeSHat()

	s.InvSHat = 1.0 / math.Max(s.SHat, 10e-10)
	s.computeFw()
	s.computeFv1()
	s.computeProduction()
	s.computeDestruction()
	s.computeCrossProduction()
	s.SourceTerm = s.Production - s.Destruction + s.CrossProduction
	return s.SourceTerm
}
