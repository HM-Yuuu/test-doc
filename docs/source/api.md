<a id="quanestimation"></a>

# quanestimation

Top-level package for quanestimation.

<a id="quanestimation.Adaptive.adaptive"></a>

# quanestimation.Adaptive.adaptive

<a id="quanestimation.Adaptive.adaptMZI"></a>

# quanestimation.Adaptive.adaptMZI

<a id="quanestimation.Adaptive"></a>

# quanestimation.Adaptive

<a id="quanestimation.AsymptoticBound.CramerRao"></a>

# quanestimation.AsymptoticBound.CramerRao

<a id="quanestimation.AsymptoticBound.CramerRao.CFIM"></a>

#### CFIM

```python
def CFIM(rho, drho, M=[], eps=1e-8)
```

Description: Calculation classical Fisher information matrix (CFIM)
             for a density matrix.

---------
Inputs
---------
rho:
    --description: parameterized density matrix.
    --type: matrix

drho:
    --description: derivatives of density matrix (rho) on all parameters
                   to be estimated. For example, drho[0] is the derivative
                   vector on the first parameter.
    --type: list (of matrix)

M:
   --description: a set of POVM. It takes the form [M1, M2, ...].
   --type: list (of matrix)

----------
Returns
----------
CFIM:
    --description: classical Fisher information matrix. If the length
                   of drho is one, the output is a float number (CFI),
                   otherwise it returns a matrix (CFIM).
    --type: float number (CFI) or matrix (CFIM)

<a id="quanestimation.AsymptoticBound.CramerRao.SLD"></a>

#### SLD

```python
def SLD(rho, drho, rep="original", eps=1e-8)
```

Description: calculation of the symmetric logarithmic derivative (SLD)
             for a density matrix.

----------
Inputs
----------
rho:
    --description: parameterized density matrix.
    --type: matrix

drho:
    --description: derivatives of density matrix (rho) on all parameters
                   to be estimated. For example, drho[0] is the derivative
                   vector on the first parameter.
    --type: list (of matrix)

rep:
    --description: the basis for the SLDs.
                   rep=original means the basis for obtained SLDs is the
                   same with the density matrix (rho).
                   rep=eigen means the SLDs are written in the eigenspace of
                   the density matrix (rho).
    --type: string {"original", "eigen"}

----------
Returns
----------
SLD:
    --description: SLD for the density matrix (rho).
    --type: list (of matrix)

<a id="quanestimation.AsymptoticBound.CramerRao.RLD"></a>

#### RLD

```python
def RLD(rho, drho, rep="original", eps=1e-8)
```

Description: calculation of the right logarithmic derivative (RLD)
             for a density matrix.

----------
Inputs
----------
rho:
    --description: parameterized density matrix.
    --type: matrix

drho:
    --description: derivatives of density matrix (rho) on all parameters
                   to be estimated. For example, drho[0] is the derivative
                   vector on the first parameter.
    --type: list (of matrix)

rep:
    --description: the basis for the RLDs.
                   rep=original means the basis for obtained RLDs is the
                   same with the density matrix (rho).
                   rep=eigen means the RLDs are written in the eigenspace of
                   the density matrix (rho).
    --type: string {"original", "eigen"}

----------
Returns
----------
RLD:
    --description: RLD for the density matrix (rho).
    --type: list (of matrix)

<a id="quanestimation.AsymptoticBound.CramerRao.LLD"></a>

#### LLD

```python
def LLD(rho, drho, rep="original", eps=1e-8)
```

Description: Calculation of the left logarithmic derivative (LLD)
            for a density matrix.

----------
Inputs
----------
rho:
    --description: parameterized density matrix.
    --type: matrix

drho:
    --description: derivatives of density matrix (rho) on all parameters
                   to be estimated. For example, drho[0] is the derivative
                   vector on the first parameter.
    --type: list (of matrix)

rep:
    --description: the basis for the LLDs.
                   rep=original means the basis for obtained LLDs is the
                   same with the density matrix (rho).
                   rep=eigen means the LLDs are written in the eigenspace of
                   the density matrix (rho).
    --type: string {"original", "eigen"}

----------
Returns
----------
LLD:
    --description: LLD for the density matrix (rho).
    --type: list (of matrix)

<a id="quanestimation.AsymptoticBound.CramerRao.QFIM"></a>

#### QFIM

```python
def QFIM(rho, drho, LDtype="SLD", exportLD=False, eps=1e-8)
```

Description: Calculation of quantum Fisher information matrix (QFIM)
            for a density matrix.

----------
Inputs
----------
rho:
    --description: parameterized density matrix.
    --type: matrix

drho:
    --description: derivatives of density matrix (rho) on all parameters
                   to be estimated. For example, drho[0] is the derivative
                   vector on the first parameter.
    --type: list (of matrix)

LDtype:
    --description: the type of logarithmic derivatives.
    --type: string {'SLD', 'RLD', 'LLD'}

exportLD:
    --description: if True, the corresponding value of logarithmic derivatives
                   will be exported.
    --type: bool

----------
Returns
----------
QFIM:
    --description: Quantum Fisher information matrix. If the length
                   of drho is 1, the output is a float number (QFI),
                   otherwise it returns a matrix (QFIM).
    --type: float number (QFI) or matrix (QFIM)

<a id="quanestimation.AsymptoticBound.CramerRao.QFIM_Bloch"></a>

#### QFIM\_Bloch

```python
def QFIM_Bloch(r, dr, eps=1e-8)
```

Description: Calculation of quantum Fisher information matrix (QFIM)
            in Bloch representation.

----------
Inputs
----------
r:
    --description: parameterized Bloch vector.
    --type: vector

dr:
    --description: derivatives of Bloch vector on all parameters to
                    be estimated. For example, dr[0] is the derivative
                    vector on the first parameter.
    --type: list (of vector)

<a id="quanestimation.AsymptoticBound.CramerRao.QFIM_Gauss"></a>

#### QFIM\_Gauss

```python
def QFIM_Gauss(R, dR, D, dD, eps=1e-8)
```

Description: Calculation of quantum Fisher information matrix (QFIM)
            for Gaussian states representation.

----------
Inputs
----------
R:

dR:

D:

dD:

<a id="quanestimation.AsymptoticBound.Holevo"></a>

# quanestimation.AsymptoticBound.Holevo

<a id="quanestimation.AsymptoticBound.Holevo.HCRB"></a>

#### HCRB

```python
def HCRB(rho, drho, W, eps=1e-8)
```

Description: solve Holevo Cramer-Rao bound via the semidefinite program (SDP).

---------
Inputs
---------
rho:
    --description: parameterized density matrix.
    --type: matrix

drho:
    --description: derivatives of density matrix (rho) on all parameters
                   to be estimated. For example, drho[0] is the derivative
                   vector on the first parameter.
    --type: list (of matrix)

W:
   --description: cost matrix.
   --type: matrix

----------
Returns
----------
Holevo bound:
    --description: the Holevo Cramer-Rao bound.
    --type: float

<a id="quanestimation.AsymptoticBound"></a>

# quanestimation.AsymptoticBound

<a id="quanestimation.BayesianBound.BayesEstimation"></a>

# quanestimation.BayesianBound.BayesEstimation

<a id="quanestimation.BayesianBound.BayesianCramerRao"></a>

# quanestimation.BayesianBound.BayesianCramerRao

<a id="quanestimation.BayesianBound.BayesianCramerRao.BCFIM"></a>

#### BCFIM

```python
def BCFIM(x, p, rho, drho, M=[], eps=1e-8)
```

Description: Calculation Bayesian version of classical Fisher information
             matrix (CFIM) for a density matrix.

---------
Inputs
---------
x:
    --description: the regimes of x for the integral.
    --type: list of arrays

p:
    --description: the prior distribution.
    --type: multidimensional array

rho:
    --description: parameterized density matrix.
    --type: multidimensional lists

drho:
    --description: derivatives of density matrix (rho) on all parameters
                   to be estimated.
    --type: multidimensional lists

M:
   --description: a set of POVM. It takes the form [M1, M2, ...].
   --type: list (of matrix)

----------
Returns
----------
BCFIM:
    --description: Bayesian version of classical Fisher information matrix.
                   If the length of x is one, the output is a float number (BCFI),
                   otherwise it returns a matrix (BCFIM).
    --type: float number (BCFI) or matrix (BCFIM)

<a id="quanestimation.BayesianBound.BayesianCramerRao.BQFIM"></a>

#### BQFIM

```python
def BQFIM(x, p, rho, drho, LDtype="SLD", eps=1e-8)
```

Description: Calculation of Bayesian version of quantum Fisher information
             matrix (QFIM) for a density matrix.

----------
Inputs
----------
x:
    --description: the regimes of x for the integral.
    --type: list of arrays

p:
    --description: the prior distribution.
    --type: multidimensional array

rho:
    --description: parameterized density matrix.
    --type: multidimensional lists

drho:
    --description: derivatives of density matrix (rho) on all parameters
                   to be estimated. For example,
    --type: multidimensional lists

LDtype:
    --description: the type of logarithmic derivatives.
    --type: string {'SLD', 'RLD', 'LLD'}

----------
Returns
----------
BQFIM:
    --description: Bayesian version of quantum Fisher information matrix.
                   If the length of x is 1, the output is a float number (BQFI),
                   otherwise it returns a matrix (BQFIM).
    --type: float number (BQFI) or matrix (BQFIM)

<a id="quanestimation.BayesianBound.ZivZakai"></a>

# quanestimation.BayesianBound.ZivZakai

<a id="quanestimation.BayesianBound"></a>

# quanestimation.BayesianBound

<a id="quanestimation.Common.common"></a>

# quanestimation.Common.common

<a id="quanestimation.Common.common.sic_povm"></a>

#### sic\_povm

```python
def sic_povm(fiducial)
```

Generate a set of POVMs by applying the $d^2$ Weyl-Heisenberg displacement operators to a
fiducial state.
The Weyl-Heisenberg displacement operators are constructioned by Fuchs et al. in the article
https://doi.org/10.3390/axioms6030021 and it is realized in QBism.

<a id="quanestimation.Common"></a>

# quanestimation.Common

<a id="quanestimation.ComprehensiveOpt.AD_Compopt"></a>

# quanestimation.ComprehensiveOpt.AD\_Compopt

<a id="quanestimation.ComprehensiveOpt.ComprehensiveStruct"></a>

# quanestimation.ComprehensiveOpt.ComprehensiveStruct

<a id="quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem"></a>

## ComprehensiveSystem

```python
class ComprehensiveSystem()
```

<a id="quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem.__init__"></a>

#### \_\_init\_\_

```python
def __init__(psi0, ctrl0, measurement0, savefile, seed, eps)
```

----------
Inputs
----------

psi0:
   --description: initial guesses of states (kets).
   --type: array

ctrl0:
    --description: initial control coefficients.
    --type: list (of vector)

measurement0:
   --description: a set of POVMs.
   --type: list (of vector)

savefile:
    --description: True: save the states (or controls, measurements) and the value of the
                         target function for each episode.
                   False: save the states (or controls, measurements) and all the value
                           of the target function for the last episode.
    --type: bool

eps:
    --description: calculation eps.
    --type: float

<a id="quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem.dynamics"></a>

#### dynamics

```python
def dynamics(tspan, H0, dH, Hc=[], ctrl=[], decay=[], ctrl_bound=[])
```

----------
Inputs
----------
tspan:
   --description: time series.
   --type: array

psi0:
   --description: initial state.
   --type: vector

measurement0:
   --description: a set of POVMs.
   --type: list (of vector)

H0:
   --description: free Hamiltonian.
   --type: matrix or a list of matrix

Hc:
   --description: control Hamiltonian.
   --type: list (of matrix)

dH:
   --description: derivatives of Hamiltonian on all parameters to
                  be estimated. For example, dH[0] is the derivative
                  vector on the first parameter.
   --type: list (of matrix)

decay:
   --description: decay operators and the corresponding decay rates.
                  decay[0][0] represent the first decay operator and
                  decay[0][1] represent the corresponding decay rate.
   --type: list

ctrl_bound:
   --description: lower and upper bounds of the control coefficients.
                  ctrl_bound[0] represent the lower bound of the control coefficients and
                  ctrl_bound[1] represent the upper bound of the control coefficients.
   --type: list

<a id="quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem.SC"></a>

#### SC

```python
def SC(W=[], M=[], target="QFIM", LDtype="SLD")
```

Description: use DE algorithm to optimize states and control coefficients.

---------
Inputs
---------
M:
    --description: a set of POVM.
    --type: list of matrix

W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem.CM"></a>

#### CM

```python
def CM(rho0, W=[])
```

Description: use DE algorithm to optimize control coefficients and the measurements.

---------
Inputs
---------
rho0:
    --description: initial state.
    --type: density matrix

W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem.SM"></a>

#### SM

```python
def SM(W=[])
```

Description: use DE algorithm to optimize states and the measurements.

---------
Inputs
---------

W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem.SCM"></a>

#### SCM

```python
def SCM(W=[])
```

Description: use DE algorithm to optimize states, control coefficients and the measurements.

---------
Inputs
---------

W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.ComprehensiveOpt.DE_Compopt"></a>

# quanestimation.ComprehensiveOpt.DE\_Compopt

<a id="quanestimation.ComprehensiveOpt.PSO_Compopt"></a>

# quanestimation.ComprehensiveOpt.PSO\_Compopt

<a id="quanestimation.ComprehensiveOpt"></a>

# quanestimation.ComprehensiveOpt

<a id="quanestimation.ControlOpt.ControlStruct"></a>

# quanestimation.ControlOpt.ControlStruct

<a id="quanestimation.ControlOpt.ControlStruct.ControlSystem"></a>

## ControlSystem

```python
class ControlSystem()
```

<a id="quanestimation.ControlOpt.ControlStruct.ControlSystem.__init__"></a>

#### \_\_init\_\_

```python
def __init__(savefile, ctrl0, load, eps)
```

----------
Inputs
----------
tspan:
   --description: time series.
   --type: array

rho0:
   --description: initial state (density matrix).
   --type: matrix

H0:
   --description: free Hamiltonian.
   --type: matrix or a list of matrix

Hc:
   --description: control Hamiltonian.
   --type: list (of matrix)

dH:
   --description: derivatives of Hamiltonian on all parameters to
                  be estimated. For example, dH[0] is the derivative
                  vector on the first parameter.
   --type: list (of matrix)

decay:
   --description: decay operators and the corresponding decay rates.
                  decay[0][0] represent the first decay operator and
                  decay[0][1] represent the corresponding decay rate.
   --type: list

ctrl_bound:
   --description: lower and upper bounds of the control coefficients.
                  ctrl_bound[0] represent the lower bound of the control coefficients and
                  ctrl_bound[1] represent the upper bound of the control coefficients.
   --type: list

savefile:
    --description: True: save the control coefficients and the value of the target function
                         for each episode.
                   False: save the control coefficients and all the value of the target
                          function for the last episode.
    --type: bool

ctrl0:
    --description: initial control coefficients.
    --type: list (of vector)

eps:
    --description: calculation eps.
    --type: float

<a id="quanestimation.ControlOpt.ControlStruct.ControlSystem.QFIM"></a>

#### QFIM

```python
def QFIM(W=[], LDtype="SLD")
```

Description: use differential evolution algorithm to update the control coefficients that maximize the
             QFI (1/Tr(WF^{-1} with F the QFIM).

---------
Inputs
---------
W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.ControlOpt.ControlStruct.ControlSystem.CFIM"></a>

#### CFIM

```python
def CFIM(M=[], W=[])
```

Description: use differential evolution algorithm to update the control coefficients that maximize the
             CFI (1/Tr(WF^{-1} with F the CFIM).

---------
Inputs
---------
M:
    --description: a set of POVM.
    --type: list of matrix

W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.ControlOpt.ControlStruct.ControlSystem.HCRB"></a>

#### HCRB

```python
def HCRB(W=[])
```

Description: use differential evolution algorithm to update the control coefficients that maximize the
             HCRB.

---------
Inputs
---------
W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.ControlOpt.DDPG_Copt"></a>

# quanestimation.ControlOpt.DDPG\_Copt

<a id="quanestimation.ControlOpt.DE_Copt"></a>

# quanestimation.ControlOpt.DE\_Copt

<a id="quanestimation.ControlOpt.GRAPE_Copt"></a>

# quanestimation.ControlOpt.GRAPE\_Copt

<a id="quanestimation.ControlOpt.PSO_Copt"></a>

# quanestimation.ControlOpt.PSO\_Copt

<a id="quanestimation.ControlOpt"></a>

# quanestimation.ControlOpt

<a id="quanestimation.Dynamics.dynamics"></a>

# quanestimation.Dynamics.dynamics

<a id="quanestimation.Dynamics.dynamics.Lindblad"></a>

## Lindblad

```python
class Lindblad()
```

General dynamics of density matrices in the form of time local Lindblad master equation.
{\partial_t \rho} = -i[H, \rho] + \sum_n {\gamma_n} {Ln.rho.Ln^{\dagger}
             -0.5(rho.Ln^{\dagger}.Ln+Ln^{\dagger}.Ln.rho)}.

<a id="quanestimation.Dynamics.dynamics.Lindblad.__init__"></a>

#### \_\_init\_\_

```python
def __init__(tspan, rho0, H0, dH, decay=[], Hc=[], ctrl=[])
```

----------
Inputs
----------
tspan:
   --description: time series.
   --type: array

rho0:
   --description: initial state (density matrix).
   --type: matrix

H0:
   --description: free Hamiltonian.
   --type: matrix

Hc:
   --description: control Hamiltonian.
   --type: list (of matrix)

dH:
   --description: derivatives of Hamiltonian on all parameters to
                  be estimated. For example, dH[0] is the derivative
                  vector on the first parameter.
   --type: list (of matrix)

ctrl:
   --description: control coefficients.
   --type: list (of array)

decay:
   --description: decay operators and the corresponding decay rates.
                  decay[0] represent a list of decay operators and
                  decay[1] represent the corresponding decay rates.
   --type: list

<a id="quanestimation.Dynamics"></a>

# quanestimation.Dynamics

<a id="quanestimation.MeasurementOpt.AD_Mopt"></a>

# quanestimation.MeasurementOpt.AD\_Mopt

<a id="quanestimation.MeasurementOpt.DE_Mopt"></a>

# quanestimation.MeasurementOpt.DE\_Mopt

<a id="quanestimation.MeasurementOpt.MeasurementStruct"></a>

# quanestimation.MeasurementOpt.MeasurementStruct

<a id="quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem"></a>

## MeasurementSystem

```python
class MeasurementSystem()
```

<a id="quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem.__init__"></a>

#### \_\_init\_\_

```python
def __init__(mtype, minput, savefile, measurement0, seed, load, eps)
```

----------
 Inputs
 ----------
savefile:
     --description: True: save the measurements and the value of the target function
                          for each episode.
                    False: save the measurements and all the value of the target
                           function for the last episode.
     --type: bool

 measurement0:
    --description: a set of POVMs.
    --type: list (of vector)

 eps:
     --description: calculation eps.
     --type: float

 notes: the Weyl-Heisenberg covariant SIC-POVM fiducial state of dimension $d$
        are download from http://www.physics.umb.edu/Research/QBism/solutions.html.

<a id="quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem.dynamics"></a>

#### dynamics

```python
def dynamics(tspan, rho0, H0, dH, Hc=[], ctrl=[], decay=[])
```

----------
Inputs
----------
tspan:
   --description: time series.
   --type: array

psi0:
    --description: initial guess of states (kets).
    --type: array

H0:
   --description: free Hamiltonian.
   --type: matrix (a list of matrix)

dH:
   --description: derivatives of Hamiltonian on all parameters to
                  be estimated. For example, dH[0] is the derivative
                  vector on the first parameter.
   --type: list (of matrix)

Hc:
   --description: control Hamiltonian.
   --type: list (of matrix)

ctrl:
    --description: control coefficients.
    --type: list (of vector)

decay:
   --description: decay operators and the corresponding decay rates.
                  decay[0][0] represent the first decay operator and
                  decay[0][1] represent the corresponding decay rate.
   --type: list

ctrl_bound:
   --description: lower and upper bounds of the control coefficients.
                  ctrl_bound[0] represent the lower bound of the control coefficients and
                  ctrl_bound[1] represent the upper bound of the control coefficients.
   --type: list

<a id="quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem.CFIM"></a>

#### CFIM

```python
def CFIM(W=[])
```

Description: use differential evolution algorithm to update the measurements that maximize the
             CFI (1/Tr(WF^{-1} with F the CFIM).

---------
Inputs
---------
W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.MeasurementOpt.PSO_Mopt"></a>

# quanestimation.MeasurementOpt.PSO\_Mopt

<a id="quanestimation.MeasurementOpt"></a>

# quanestimation.MeasurementOpt

<a id="quanestimation.Resources.Resources"></a>

# quanestimation.Resources.Resources

<a id="quanestimation.Resources.Resources.SpinSqueezing"></a>

#### SpinSqueezing

```python
def SpinSqueezing(rho, basis="Dicke", output="KU")
```

Description: Calculate spin squeezing parameter for a density matrix.

---------
Inputs
---------
N :
   --description: particle number.
   --type: int
rho :
   --description: density matrix.
   --type: matrix

output:
   --description: if output=='KU',the output of the squeezing parameter is defined by
                  Kitagawa and Ueda. if output=='WBIMH',the output of the squeezing
                  parameter is defined by Wineland et al.
   --type: string {'KU', 'WBIMH'}

----------
Returns
----------
Xi :
   --description: squeezing_parameter.
   --type: float

<a id="quanestimation.Resources"></a>

# quanestimation.Resources

<a id="quanestimation.StateOpt.AD_Sopt"></a>

# quanestimation.StateOpt.AD\_Sopt

<a id="quanestimation.StateOpt.DDPG_Sopt"></a>

# quanestimation.StateOpt.DDPG\_Sopt

<a id="quanestimation.StateOpt.DE_Sopt"></a>

# quanestimation.StateOpt.DE\_Sopt

<a id="quanestimation.StateOpt.NM_Sopt"></a>

# quanestimation.StateOpt.NM\_Sopt

<a id="quanestimation.StateOpt.PSO_Sopt"></a>

# quanestimation.StateOpt.PSO\_Sopt

<a id="quanestimation.StateOpt.StateStruct"></a>

# quanestimation.StateOpt.StateStruct

<a id="quanestimation.StateOpt.StateStruct.StateSystem"></a>

## StateSystem

```python
class StateSystem()
```

<a id="quanestimation.StateOpt.StateStruct.StateSystem.__init__"></a>

#### \_\_init\_\_

```python
def __init__(savefile, psi0, seed, load, eps)
```

----------
Inputs
----------
savefile:
    --description: True: save the states and the value of target function for each episode .
                   False: save the states and the value of target function for the last episode.
    --type: bool

psi0:
    --description: initial guess of states (kets).
    --type: array

seed:
    --description: random seed.
    --type: int

eps:
    --description: machine eps.
    --type: float

<a id="quanestimation.StateOpt.StateStruct.StateSystem.dynamics"></a>

#### dynamics

```python
def dynamics(tspan, H0, dH, Hc=[], ctrl=[], decay=[])
```

----------
Inputs
----------
tspan:
   --description: time series.
   --type: array

H0:
   --description: free Hamiltonian.
   --type: matrix (a list of matrix)

dH:
   --description: derivatives of Hamiltonian on all parameters to
                  be estimated. For example, dH[0] is the derivative
                  vector on the first parameter.
   --type: list (of matrix)

Hc:
   --description: control Hamiltonian.
   --type: list (of matrix)

ctrl:
    --description: control coefficients.
    --type: list (of vector)

decay:
   --description: decay operators and the corresponding decay rates.
                  decay[0][0] represent the first decay operator and
                  decay[0][1] represent the corresponding decay rate.
   --type: list

ctrl_bound:
   --description: lower and upper bounds of the control coefficients.
                  ctrl_bound[0] represent the lower bound of the control coefficients and
                  ctrl_bound[1] represent the upper bound of the control coefficients.
   --type: list

<a id="quanestimation.StateOpt.StateStruct.StateSystem.QFIM"></a>

#### QFIM

```python
def QFIM(W=[], LDtype="SLD")
```

Description: use autodifferential algorithm to search the optimal initial state that maximize the
             QFI (1/Tr(WF^{-1} with F the QFIM).

---------
Inputs
---------
W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.StateOpt.StateStruct.StateSystem.CFIM"></a>

#### CFIM

```python
def CFIM(M=[], W=[])
```

Description: use autodifferential algorithm to search the optimal initial state that maximize the
             CFI (1/Tr(WF^{-1} with F the CFIM).

---------
Inputs
---------
M:
    --description: a set of POVM.
    --type: list of matrix

W:
    --description: weight matrix.
    --type: matrix

<a id="quanestimation.StateOpt"></a>

# quanestimation.StateOpt

