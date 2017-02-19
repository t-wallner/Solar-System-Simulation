/*
 * 
 * Computational Methods, Mini-Project: Particle3D class
 *
 * Themis Rallis, Tom Emmens, Thomas Wallner 2014
 *
 *
 */

// Imports
import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;
import java.lang.Math;

public class Particle3D
{

	// Variables
	private double mass; // mass of particle
	private String name; // name of particle
	Vector3D position = new Vector3D(); // position vector of particle
	Vector3D velocity = new Vector3D(); // velocity vector of particle

	// Constructors
	// Default Constructor
	public Particle3D()
	{
		mass = 1.0;
		name = "Nameless";
		velocity.setVector3D(0, 0, 0);
		position.setVector3D(0, 0, 0);
	}

	// Explicit Constructor
	public Particle3D(String name, double mass, Vector3D p, Vector3D v)
	{
		this.name = name;
		this.mass = mass;
		this.position = p;
		this.velocity = v;
	}

	// Methods
	// Instance Methods

	// Return the name and position of a particle in a VMD friendly format
	public String toString()
	{
		return (name + " " + position.getX() + " " + position.getY() + " " + position.getZ());
	}

	// Returns a particle with characteristics extracted from a given input file
	public Particle3D scan(String filename) throws IOException
	{
		// Read input from file
		BufferedReader file = new BufferedReader(new FileReader(filename));
		Scanner scan = new Scanner(file);

		// Allocate tokens
		name = scan.next();
		mass = scan.nextDouble();
		position.setVector3D(scan.nextDouble(), scan.nextDouble(), scan.nextDouble());
		velocity.setVector3D(scan.nextDouble(), scan.nextDouble(), scan.nextDouble());
		scan.close();
		return new Particle3D(name, mass, position, velocity);
	}

	// Returns the Kintetic energy
	public double kineticEnergy()
	{
		return 0.5 * mass * velocity.magnitudeSq();
	}

	// Integration support - evolve the position : dx = v * dt
	public void leapPosition(double dt)
	{
		position.setVector3D(position.getX() + velocity.getX() * dt,
				position.getY() + velocity.getY() * dt, position.getZ() + velocity.getZ() * dt);
	}

	// Integration support - evolve the velocity : dv = f/m * dt
	public void leapVelocity(double dt, Vector3D force)
	{
		velocity.setVector3D(velocity.getX() + force.getX() * dt / mass,
				velocity.getY() + force.getY() * dt / mass, velocity.getZ() + force.getZ() * dt
						/ mass);
	}

	// Getters and Setters
	public String getName()
	{
		return name;
	}

	public void setName(String newname)
	{
		name = newname;
	}

	public double getMass()
	{
		return mass;
	}

	public void setMass(double m)
	{
		this.mass = m;
	}

	public Vector3D getPosition()
	{
		return position;
	}

	public Vector3D getVelocity()
	{
		return velocity;
	}

	public void setVelocity(Vector3D v)
	{
		this.velocity = new Vector3D(v);
	}

	// Static Methods
	
	// Returns the magntitude of the separation between two particles
	public static double separation(Particle3D a, Particle3D b)
	{
		return Vector3D.subtract(a.getPosition(), b.getPosition()).magnitude();
	}

	// Integration support - evolve all positions : dx = v * dt
	public static void leapPositionAll(ArrayList<Particle3D> particles, double timestep)
	{
		for (Particle3D p : particles)
		{
			p.leapPosition(timestep);
		}
	}

	// Integration support - evolve all velocities : dv = f/m * dt
	public static void leapVelocityAll(ArrayList<Particle3D> particles, double timestep)
	{
		for (Particle3D p : particles)
		{
			p.leapVelocity(timestep, new Vector3D(resultantForce(particles, p)));
		}
	}

	// Returns the force vector between two particles
	public static Vector3D forceBetweenParticles(Particle3D a, Particle3D b)
	{
		double G = 8.89 * Math.pow(10, -10); // Gravitational constant in terms of AU and Earth mass.
		return Vector3D.subtract(b.getPosition(), a.getPosition()).scale(
				(G * a.getMass() * b.getMass() / Math.pow(Particle3D.separation(a, b), 3)));
	}

	// Returns the resultant force vector on the focusParticle due to all the others in the system
	public static Vector3D resultantForce(ArrayList<Particle3D> particles, Particle3D focusParticle)
	{
		ArrayList<Particle3D> systemCopy = new ArrayList<Particle3D>(particles);
		Vector3D force = new Vector3D();
		systemCopy.remove(focusParticle);
		for (Particle3D p : systemCopy)
		{
			force = new Vector3D(Vector3D.add(forceBetweenParticles(focusParticle, p), force));
		}
		return force;
	}

	// Returns the total system kinetic energy
	public static double systemKineticEnergy(ArrayList<Particle3D> particles)
	{
		double totalKineticEnergy = 0;
		for (Particle3D p : particles)
		{
			totalKineticEnergy += p.kineticEnergy();
		}
		return totalKineticEnergy;
	}

	// Returns the potential energy between two particles
	public static double potentialEnergyBetweenParticles(Particle3D a, Particle3D b)
	{
		return (-a.getMass() * b.getMass() / Particle3D.separation(a, b));
	}

	//  Returns the total system potential energy
	public static double systemPotentialEnergy(ArrayList<Particle3D> particles)
	{
		double systemPotential = 0;
		for (int i = 0; i < particles.size() - 1; i++)
		{
			for (int j = i + 1; j < particles.size(); j++)
			{
				systemPotential += potentialEnergyBetweenParticles(particles.get(i),particles.get(j));
			}
		}
		return systemPotential;
	}

	// Returns the total system energy : PE + KE
	public static double systemTotalEnergy(ArrayList<Particle3D> particles)
	{
		return systemPotentialEnergy(particles) + systemKineticEnergy(particles);
	}

	// Returns an array list of particles with characteristics extracted from a given input file
	public static ArrayList<Particle3D> scanAllParticles(String filename) throws IOException
	{
		ArrayList<Particle3D> particles = new ArrayList<Particle3D>();
		// Read input from file
		BufferedReader file = new BufferedReader(new FileReader(filename));
		Scanner scan = new Scanner(file);

		while (scan.hasNext())
		{
			// Allocate tokens
			String name = scan.next();
			double mass = scan.nextDouble();
			Vector3D position = new Vector3D(scan.nextDouble(), scan.nextDouble(),
					scan.nextDouble());
			Vector3D velocity = new Vector3D(scan.nextDouble(), scan.nextDouble(),
					scan.nextDouble());
			particles.add(new Particle3D(name, mass, position, velocity));
		}

		scan.close();

		return particles;
	}

	// Calculates the initial velocity of the Sun such that the total system momentum is zero
	public static void zeroSystemMomentum(ArrayList<Particle3D> particles, Particle3D focusParticle)
	{
		ArrayList<Particle3D> systemCopy = new ArrayList<Particle3D>(particles);
		systemCopy.remove(focusParticle);
		Vector3D momentum = new Vector3D();
		for (Particle3D p : systemCopy)
		{
			momentum = Vector3D.add(p.getVelocity().scale(p.getMass()), momentum);
		}

		Vector3D adjustedVelocity = momentum.scale(-1.0 / focusParticle.getMass());

		focusParticle.setVelocity(adjustedVelocity);
	}

	// Calculates the change in angle of a particle in the x-y plane, and then updates the orbit
	public static void updateTheta(ArrayList<Particle3D> particles, double[] xposition,
					double[] yposition, double[] changeTheta, int[] orbits)
	{
		double angle1;
		double angle2;
		for (int i = 0; i < particles.size(); i++)
		{
			// update period for all planets
			if (!particles.get(i).getName().equals("Moon"))
			{
				// calculate angle based on previous position
				angle1 = Math.abs(Math.atan2(yposition[i], xposition[i]));
				// calculate angle based on current position
				angle2 = Math.abs(Math.atan2(particles.get(i).getPosition().getY(), particles
						.get(i).getPosition().getX()));
				// find change in angle (theta)
				changeTheta[i] = changeTheta[i] + Math.abs(angle2 - angle1);
				// only increase orbit once the total change in angle (theta)
				// exceeds 2pie
				if (changeTheta[i] > 2 * Math.PI)
				{
					orbits[i] += 1;
					changeTheta[i] = changeTheta[i] - 2 * Math.PI;
				}
			}
			// update period for moon
			else
			{
				// calculate angle based on previous position
				angle1 = Math.abs(Math.atan2(Math.abs(yposition[i] - yposition[i - 1]),
						Math.abs(xposition[i] - xposition[i - 1])));
				// calculate angle based on current position
				angle2 = Math.abs(Math.atan2(
						Math.abs(particles.get(i).getPosition().getY()
								- particles.get(i - 1).getPosition().getY()),
						Math.abs(particles.get(i).getPosition().getX()
								- particles.get(i - 1).getPosition().getX())));
				// find change in angle (theta)
				changeTheta[i] = changeTheta[i] + Math.abs(angle2 - angle1);
				// only increase orbit once the total change in angle (theta)
				// exceeds 2pie
				if (changeTheta[i] > 2 * Math.PI)
				{
					orbits[i] += 1;
					changeTheta[i] = changeTheta[i] - 2 * Math.PI;
				}
			}
		}
	}

	// Updates the closest point of a particle to the sun
	public static void updateClosestPoints(ArrayList<Particle3D> particles, double[] closestPoints, int sunIndex)
	{
		for (int i = 0; i < closestPoints.length; i++)
		{
			if (Particle3D.separation(particles.get(sunIndex), particles.get(i)) < closestPoints[i])
			{
				closestPoints[i] = Particle3D.separation(particles.get(sunIndex), particles.get(i));
			}
		}
	}	
	
	// Updates the farthest point of a particle from the sun 
	public static void updateFarthestPoints(ArrayList<Particle3D> particles, double[] farthestPoints, int sunIndex)
	{
		for (int i = 0; i < farthestPoints.length; i++)
		{	
			if (Particle3D.separation(particles.get(sunIndex), particles.get(i)) > farthestPoints[i])
			{ // farthest point from sun
				farthestPoints[i] = Particle3D.separation(particles.get(sunIndex), particles.get(i));
			}

		}
	}
}
