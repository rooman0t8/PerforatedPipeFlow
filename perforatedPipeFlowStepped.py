#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Andrew Hedrick
Description:
   Solve flow parameters for a perforated pipe.
   See the inline comments for details.
   All parameters are defined in source.
"""

import math

# Problem statement
# Air flows through a PVC pipe with a given diameter and length which is capped at the end. The pipe
# is perforated with holes of equal diameter at an even spacing. Pressure outside the
# pipe is sea-level standard air.  A fan entrance at point 1 draws air from the
# pipe with known performance.  Estimate the pressure at point 2 and if there is any
# corresponding change in the flow rate at the holes.
#
#   ---------------------------------------------------|
# <---                                                 | d
#   ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -|
#   1                 L                                2

# Assumptions
# Flow is incompressible and steady state.
# Neglecting minor pipe losses for initial solution.
# Due to endcap velocity at point 2 is zero.  The total inflow from the holes must equal 
#   flow rate at point 1 to satisfy conservation of mass.
# Holes are approximated as thin plate orifices.

# References
# Ref1, Fluid Mechanics Fith Ed. by Frank M. White, 2003
# Specific urls given in-line

# Programing notes
# I have deliberately avoided importing UnitsNet or other packages to minimize dependencies for sharing.
# All calculations in SI units 



def calculateFlowStep(velocityA, pressureA, lengthStep, areaPerStep, pipe):
   """
   Calculates the flow through a control volume account for friction losses and perforation inflow.
   """
   # Calculate velocity and Reynolds number
   # ---------------------------------
   flowRateA = velocityA * pipe.area
   #print("flowRateA = ", flowRateA, "m^3/s")

   reNum = (density * velocityA * pipe.diameter) / kViscosity
   #print("Renolds number = ", reNum)

   flowRateHoles = Cd * areaPerStep * (2 * math.fabs(pressureAtm - pressureA)) ** 0.5
   #print("flow rate holes =", flowRateHoles, "m^3/s")

   flowRateB = flowRateA - flowRateHoles
   velocityB = flowRateB / pipe.area

   # Calculate friction components
   # ---------------------------------
   # From multiple online resources setting the absolute roughness coefficient for PVC
   # https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html
   # https://www.ti-soft.com/en/support/help/heatingdesign/project/parameters/hydronic-settings/pipe-roughness/absolute_pipe_roughness
   # https://www.pipeflow.com/pipe-pressure-drop-calculations/pipe-roughness
   roughness = 0.0015 # mm

   # Solving the Swamee equation for the Darcy–Weisbach friction factor; it is selected for
   # validity across the entire flow regime as the Renolds number drops.
   # https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Swamee_equation
   A = (64/reNum) ** 8
   B = math.log((roughness / (3.7 * pipe.diameter)) + (5.74 / (reNum ** 0.9)))
   C = (2500 / reNum) ** 6
   frictionFactor = (A + 9.5 * ((B - C) ** -16)) ** (1/8)
   #print("frictionFactor =", frictionFactor)

   # Solving for the wall shear stress using Eq 6.11 from Ref1
   shearStressWall = (frictionFactor / 8) * density * (velocityA ** 2)
   #print("shearStressWall = ", shearStressWall, "Pascals")


   # Solving Control Volume 
   # ---------------------------------
   # Assigning a control volume to enclose points 1, 2 and the pipe walls.
   # The linear momentum equation based on Eq 3.40 from Ref1
   # Note the positive term for v1 as an outlet and positive term for friction
   # to oppose the flow direction.
   # The original equation for readability
   #   sumForces = p1A - p2A + tau*Pi*D*L = -rho*A*V2^2 + rho*A*V1^2

   # deltaP = pA - pB
   deltaP = -shearStressWall * ((math.pi * pipe.diameter * lengthStep) / pipe.area) + density * velocityA ** 2 - density * velocityB ** 2
   #print("pA - pB = ", deltaP, "Pascals")

   pressureB = pressureA - deltaP
   #print("pressureB = ", pressureB, "Pascals")
   return velocityB, pressureB


def solution(flowState, areaHolesTotal, pipe):
   """
   Estimates flow through the pipe for a given total area for the holes. 
   """

   velocityA = flowState.velocity1
   pressureA = flowState.pressure1

   # Steps is set to balance grainularity with a size to have an actual effect.
   # No care has been given to optimize or safeguard the solution against very
   # small values used in calculations.
   steps = 500
   lastStep = steps - 1
   
   # The target velocity for the last step to approximate the known condition
   # velocity2 = 0.  Given the steps and area increments this value should be adjusted;
   # too small a value can prevent converging on a solution.
   targetVelocity = 0.1
   areaPerStep = areaHolesTotal / steps
   lengthStep = pipe.length / steps

   flowState.starved = False
   for step in range(steps):
      velocityB, pressureB = calculateFlowStep(velocityA, pressureA, lengthStep, areaPerStep, pipe)
      if velocityB < 0:
            flowState.starved = True
            #print("! Early dead flow velocityB =", velocityB, "at step ", step)
            break
      if velocityB < targetVelocity:
         if step == lastStep:
            flowState.solved = True

      velocityA = velocityB
      pressureA = pressureB

   flowState.velocity2 = velocityB
   #print("flowState.velocity2 = ", flowState.velocity2)
   flowState.pressure2 = pressureB
   #print("flowState.pressure2 = ", flowState.pressure2)

   return flowState


class flowState():
   velocity1 = 0
   pressure1 = 0
   velocity2 = 0
   pressure2 = 0
   solved = False
   starved = False
   steps = []

class pipe():
   diameter = 0
   area= 0
   length = 0


def main():
   # Known Parameters
   # ---------------------------------
   # Fan model Rn4EC-4
   # https://www.fantech.net/en-us/products/fans-and-accessories/radon-fans/rn?sku=99923
   flowRate1Cfm = 175
   flowRate1 = flowRate1Cfm * cfmToCubicMetersPerSecond
   print("flowRate1 =", flowRate1, "(", flowRate1Cfm, ")", "m^3/s (cfm)")

   pressureFanInWC = -2.25
   pressureFan = pressureFanInWC * inchesWaterColumnToPascals
   print("pressureFan =", pressureFan, "(", pressureFanInWC, ")", "Pascals (InWC)")

   diameterPipeIn = 4
   pipe.diameter = diameterPipeIn * inchsToMeters
   pipe.area = (math.pi / 4) * pipe.diameter ** 2

   lengthPipeFt = 100
   pipe.length = lengthPipeFt * feetToMeters
   print("lengthPipe =", pipe.length, "meters")

   flowState.velocity1 = flowRate1 / pipe.area
   flowState.pressure1 = pressureAtm + pressureFan
   print("pressure1 = ", flowState.pressure1, "Pascals")

   # An initial estimate of the hole area to begin iterating.
   # A first guess can be calculated below.
      # Using Eq 6.104 Ref1
      # Assuming beta = 0 for small holes
      # deltaP = 100 # first pass estimate
      #areaHolesTotal = flowRate1 / (Cd * ((2 * deltaP / density) ** 0.5))
   areaHolesTotal = 0.005
   
   # The change in the area should be sized no smaller than a reasonable drill bit
   increment = 0.00001

   # Start Solution
   while flowState.solved == False:
      if flowState.starved:
         areaHolesTotal -= increment
      else:
         areaHolesTotal += increment
      #print("areaHolesTotal = ", areaHolesTotal, "m^2")

      solution(flowState, areaHolesTotal, pipe)
      #print("flowState.solved =", flowState.solved)

   print()
   print("-------------------------")
   print("areaHolesTotal = ", areaHolesTotal, "m^2")
   print("velocity1 = ", flowState.velocity1, "m/s")
   print("pressure1 = ", flowState.pressure1, "Pascals")

   print("velocity2 = ", flowState.velocity2, "m/s")
   print("pressure2 =", flowState.pressure2, "Pascals")

   deltaP = flowState.pressure2 - flowState.pressure1
   print("  deltaP = ", deltaP, " percent of fan pressure %", 100 *(deltaP / pressureFan))


# ---------------------------------
# Conversion Constants
feetToMeters = 0.3048
inchsToMeters = 0.0254
cfmToCubicMetersPerSecond = 0.0004719
inchesWaterColumnToPascals = 248.84

# Standard air at sea-level
pressureAtm = 101350 # Pascals
density = 1.225 # kg/m^3
kViscosity = 1.789E-5 # Pa*s
Cd = 0.6 # Assumption based on examples in Ref1

if __name__ == '__main__':
    main()
