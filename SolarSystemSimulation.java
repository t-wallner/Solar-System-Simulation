/*
 * 
 * Computational Methods, Mini-Project: SolarSystemSimulation class (main)
 *
 * Themis Rallis, Tom Emmens, Thomas Wallner 2014
 *
 *
 */

// Imports
import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;

public class SolarSystemSimulation
{

	public static void main(String[] args) throws IOException // IO Exception to catch Input/Output Errors
	{

		// Create array of particles and scan particle characteristics from file
		ArrayList<Particle3D> particles = Particle3D.scanAllParticles(args[0]);

		// Initial & final energy
		double initialEnergy = Particle3D.systemTotalEnergy(particles);
		double finalEnergy = Particle3D.systemTotalEnergy(particles);
		// Initial time
		double time = 0;
		// Array to hold number of orbits around the sun
		int[] orbits = new int[particles.size()];
		// Array to hold the realistic year lengths of planets
		double[] yearLengths = new double[particles.size()];

		// Input time step length and number of time steps from file
		BufferedReader file = new BufferedReader(new FileReader(args[1]));
		Scanner scan = new Scanner(file);
		double timestep = scan.nextDouble();
		int steps = scan.nextInt();
		scan.close();

		// Open an output file for trajectories.xyz
		PrintWriter outputTrajectories = new PrintWriter(new FileWriter(args[2]));
		// Open an output file for results (i.e. Orbit counts, year lengths,
		// aphelions, perihelions, initial to final energy ratio)
		PrintWriter outputResults = new PrintWriter(new FileWriter(args[3]));

		// Set the velocity of the sun so that the solar system has zero total
		// momentum
		Particle3D.zeroSystemMomentum(particles, particles.get(0));

		// Print the initial position of each particle to a file
		outputTrajectories.println(particles.size());
		outputTrajectories.println("Point = 1");
		for (Particle3D p : particles)
		{
			outputTrajectories.println(p.toString());
		}

		// Leap velocity with half time step to get velocity at (t + dt/2) as
		// required by by the leap frog algorithm
		Particle3D.leapVelocityAll(particles, timestep);

		// Initiate arrays to store the closest and farthest distance each planet reaches from the sun.
		double[] closestPoints = new double[particles.size()]; // used for perihelion if particle completes a full orbit
		double[] farthestPoints = new double[particles.size()]; // used for aphelion if particle completes a full orbit
		

		// Store the initial separation between each particle and the sun
		for (int i = 0; i < particles.size(); i++)
		{
			closestPoints[i] = Particle3D.separation(particles.get(0), particles.get(i));
			farthestPoints[i] = Particle3D.separation(particles.get(0), particles.get(i));
		}

		// Create array for velocities to track number of orbits
		double[] xPositions = new double[particles.size()];
		double[] yPositions = new double[particles.size()];
		double[] changeInTheta = new double[particles.size()];

		// Start leap frog loop
		for (int i = 0; i < steps; i++)
		{
			// Increase time
			time += timestep;

			// Store positions of all particles before the position leap.
			for (int j = 0; j < particles.size(); j++)
			{
				xPositions[j] = particles.get(j).getPosition().getX();
				yPositions[j] = particles.get(j).getPosition().getY();
			}

			// Update position.
			Particle3D.leapPositionAll(particles, timestep);

			// Update closestPoints and farthestPoints arrays.			
			Particle3D.updateFarthestPoints(particles, farthestPoints, 0);
			Particle3D.updateClosestPoints(particles, closestPoints, 0);

			// Update velocity
			// (The leapVelocityAll method calculates the particle specific force, then leaps the velocity)
			Particle3D.leapVelocityAll(particles, timestep);

			// Update the number of full orbits around the sun
			Particle3D.updateTheta(particles, xPositions, yPositions, changeInTheta, orbits);

			// Record realistic year length of first using complete orbit.
			for (int k = 0; k < particles.size(); k++)
			{
				// Year length recorded at first step after the particle has completed a full orbit.
				// Only assigned once (On the condition it is currently unassigned). 
				if (orbits[k] == 1 && yearLengths[k] == 0) { yearLengths[k] = time; };
				// Particles that do not complete a full orbit will have the default year Length -1.0
			}

			// Print results to file in VMD friendly format
			outputTrajectories.println(particles.size());
			outputTrajectories.println("Point = " + (2 + i));
			for (int l = 0; l < particles.size(); l++)
			{
				outputTrajectories.println(particles.get(l).toString());
			}

		}
		
		
		// Initiate two arrays to store aphelion and perihelion distances for each particle.
		double[] aphelions = new double[particles.size()];
		double[] perihelions = new double[particles.size()];
		
		// If a particle has not completed a full orbit during simulation assign its aphelion and perihelion values to -1.0
		// otherwise assign the closest and farthest points from the sun.
		for (int i = 0; i < particles.size(); i++)
		{	
			aphelions[i] = (orbits[i] > 0) ? farthestPoints[i] : -1.0;
			perihelions[i] = (orbits[i] > 0) ? closestPoints[i] : -1.0;
		}
		
		// Print results to results file
		outputResults.println("Number of full orbits: ");
		for (int i = 0; i < particles.size(); i++)
		{
			outputResults.println(particles.get(i).getName() + " " + orbits[i]);
		}
		// Print out the aphelion and perihelion
		outputResults.println("\nCalculated Aphelions & Perihelions:");
		for (int i = 0; i < particles.size(); i++)
		{
			outputResults.println(particles.get(i).getName() + " Aphelion: " + aphelions[i]
					+ " Perihelion: " + perihelions[i]);
		}
		// Print out the year lengths
		outputResults.println("\nAccurate Year Length (in days):");
		for (int i= 0; i< particles.size(); i++)
		{
			if (yearLengths[i] != 0)
			{
				outputResults.println(particles.get(i).getName() + " " + yearLengths[i]);

			} else
			{
				outputResults.println(particles.get(i).getName()+" has not completed an orbit!");
			}
			
		}
		// Print out the year ratios
		outputResults.println("\nYear ratios:");
		// only proceed if earth is part of simulation
		if(particles.size()>3)
		{
			for (int i = 0; i < particles.size(); i++)
			{
				// Calculate year ratio with respect to the third input planet
				// (MUST BE THE EARTH!)
				if (yearLengths[3] != 0)
				{
					outputResults.println(particles.get(i).getName() + " " + yearLengths[i]
							/ yearLengths[3]);
				} else
				{
					outputResults.println("The earth has not completed an orbit!");
					break;
				}
			}
		}
		// Print out the ratio of final and initial energies
		outputResults.println("\nRatio of final and initial system energy:");
		
		// Update final Energy		
		finalEnergy = Particle3D.systemTotalEnergy(particles);
		// Print out the ratio final Energy to initial Energy
		outputResults.println(finalEnergy / initialEnergy);
		
		// Close writer files
		outputTrajectories.close();
		outputResults.close();

	}
}
