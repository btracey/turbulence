package sst

import (
	"math"

	"github.com/btracey/fluid"
	"github.com/btracey/fluid/fluid2d"
	"github.com/gonum/unit"
)

const (
	SigmaK1      = 0.85
	SigmaK2      = 1.0
	BetaStar     = 0.09
	SqrtBetaStar = 0.3
	SigmaW1      = 0.5
	SigmaW2      = 0.856
	Kappa        = 0.41
	Beta1        = 0.075
	Beta2        = 0.0828
	A1           = 0.31
	Gamma1       = Beta1/BetaStar - SigmaW1*Kappa*Kappa/SqrtBetaStar
	Gamma2       = Beta2/BetaStar - SigmaW2*Kappa*Kappa/SqrtBetaStar
)

type SST2Cache struct {
	// Inputs
	VelGrad   fluid2d.VelGrad // Velocity gradient tensor du_i / dx_j
	WallDist  unit.Length
	K         fluid.KineticEnergy
	Nu        fluid.KinematicViscosity
	Rho       fluid.Density
	OmegaDiss fluid.SpecificDissipation
	DKDX      [2]float64
	DOmegaDX  [2]float64

	// Cached
	S                    fluid2d.StressTensor
	OmegaVort            fluid.Vorticity
	Arg2                 float64
	F2                   float64
	MuT                  fluid.Viscosity
	Tau                  fluid2d.ReynoldsStress
	P                    float64
	CDkw                 float64
	Arg1                 float64
	F1                   float64
	ProductionOmega      float64
	DestructionOmega     float64
	CrossProductionOmega float64

	// Final results
	SourceK     float64
	SourceOmega float64

	// Nondimensionalized
	NondimSourceK     float64
	NondimSourceOmega float64
}

func (s *SST2Cache) Compute() {
	s.S = s.VelGrad.SymmetricPart()
	//asym :=
	s.OmegaVort = s.VelGrad.Vorticity()
	s.Arg2 = Arg2(s.K, s.WallDist, s.Nu, s.OmegaDiss)
	s.F2 = F2(s.Arg2)

	s.MuT = TurbVisc(s.Rho, s.K, s.OmegaDiss, s.OmegaVort, s.F2)
	s.Tau = ReStress(s.MuT, s.S, s.K, s.Rho)
	s.P = Production(s.Tau, s.VelGrad)
	s.CDkw = CDkw(s.Rho, s.OmegaDiss, s.DKDX, s.DOmegaDX)
	s.Arg1 = Arg1(s.K, s.OmegaDiss, s.WallDist, s.Nu, s.Rho, s.CDkw)
	s.F1 = F1(s.Arg1)
	s.ProductionOmega = ProductionOmega(s.F1, s.Rho, s.MuT, s.P)
	s.DestructionOmega = DestructionOmega(s.F1, s.Rho, s.OmegaDiss)
	s.CrossProductionOmega = CrossProductionOmega(s.F1, s.Rho, s.OmegaDiss, s.DKDX, s.DOmegaDX)

	s.SourceK = SourceK(s.P, s.Rho, s.OmegaDiss, s.K)
	s.SourceOmega = SourceOmega(s.ProductionOmega, s.DestructionOmega, s.CrossProductionOmega)

	// Don't want to divide directly because of numerical errors
	kNondim := float64(s.Rho) * float64(s.K) * float64(s.OmegaDiss)
	s.NondimSourceK = s.P/kNondim - BetaStar
	omegaNondim := float64(s.Rho) * float64(s.OmegaDiss) * float64(s.OmegaDiss)
	beta := Blend(s.F1, Beta1, Beta2)
	s.NondimSourceOmega = s.ProductionOmega/omegaNondim - beta + s.CrossProductionOmega/omegaNondim
}

// Can be used for
type SST2 struct {
}

func F1(arg1 float64) float64 {
	return math.Tanh(arg1 * arg1 * arg1 * arg1)
}

func F2(arg2 float64) float64 {
	return math.Tanh(arg2 * arg2)
}

func Arg2(k fluid.KineticEnergy, d unit.Length, nu fluid.KinematicViscosity, omega fluid.SpecificDissipation) float64 {
	t1 := 2 * math.Sqrt(float64(k)) / (BetaStar * float64(omega) * float64(d))
	t2 := 500 * float64(nu) / (float64(d) * float64(d) * float64(omega))
	return math.Max(t1, t2)
}

func Production(tau fluid2d.ReynoldsStress, u fluid2d.VelGrad) float64 {
	return tau.UU()*u.DUDX() + tau.UV()*u.DUDY() + tau.UV()*u.DVDX() + tau.VV()*u.DVDY()
}

func ReStress(mut fluid.Viscosity, S fluid2d.StressTensor, k fluid.KineticEnergy, rho fluid.Density) fluid2d.ReynoldsStress {
	// equation says du_k/dx_k, but the trace of s is the same.
	trace := S.At(0, 0) + S.At(1, 1)
	tau := fluid2d.ReynoldsStress{}
	for i := 0; i < 2; i++ {
		for j := i; j < 2; j++ {
			term := float64(mut) * 2 * S.At(i, j)
			if i == j {
				term -= float64(mut) * 2.0 / 3.0 * trace
				term -= 2.0 / 3.0 * float64(rho) * float64(k)
				tau.Set(i, j, term)
			}
		}
	}
	return tau
}

func TurbVisc(rho fluid.Density, k fluid.KineticEnergy, omega fluid.SpecificDissipation, Omega fluid.Vorticity, F2 float64) fluid.Viscosity {
	return fluid.Viscosity(float64(rho) * A1 * float64(k) / math.Max(A1*float64(omega), F2*float64(Omega)))
}

func SourceK(P float64, rho fluid.Density, omega fluid.SpecificDissipation, k fluid.KineticEnergy) float64 {
	return P - BetaStar*float64(rho)*float64(omega)*float64(k)
}

func CDkw(rho fluid.Density, omega fluid.SpecificDissipation, dkdx, dOmegaDx [2]float64) float64 {
	cp := dkdx[0]*dOmegaDx[0] + dkdx[1]*dOmegaDx[1]
	term := 2 * float64(rho) * SigmaW2 * cp / float64(omega)
	return math.Max(term, 1e-20)
}

func Arg1(k fluid.KineticEnergy, omega fluid.SpecificDissipation, d unit.Length,
	nu fluid.KinematicViscosity, rho fluid.Density, cdkw float64) float64 {
	d2 := float64(d) * float64(d)
	term1 := math.Sqrt(float64(k)) / (BetaStar * float64(omega) * float64(d))
	term2 := 500 * float64(nu) / (d2 * float64(omega))
	term3 := 4 * float64(rho) * SigmaW2 * float64(k) / (cdkw * d2)
	return math.Min(math.Max(term1, term2), term3)
}

func Blend(f1, one, two float64) float64 {
	return f1*one + (1-f1)*two
}

func ProductionOmega(f1 float64, rho fluid.Density, muT fluid.Viscosity, P float64) float64 {
	gamma := Blend(f1, Gamma1, Gamma2)
	nuT := float64(muT) / float64(rho)
	return (gamma / nuT) * P
}

func DestructionOmega(f1 float64, rho fluid.Density, omega fluid.SpecificDissipation) float64 {
	beta := Blend(f1, Beta1, Beta2)
	return beta * float64(rho) * float64(omega) * float64(omega)
}

func CrossProductionOmega(f1 float64, rho fluid.Density, omega fluid.SpecificDissipation, dkdx [2]float64, dOmegaDX [2]float64) float64 {
	prod := dkdx[0]*dOmegaDX[0] + dkdx[1]*dOmegaDX[1]
	return 2 * (1 - f1) * float64(rho) * SigmaW2 * prod / float64(omega)
}

func SourceOmega(prod, dest, crossProd float64) float64 {
	return prod - dest + crossProd
}
