/*
 * 
 * Computational Methods, Exercise 2: Vector 3D class
 *
 * Thomas Wallner, 2014
 */ 

public class Vector3D {

    // Variables
    private double x,y,z; // Coordinates of Vector

    // Constructors
    // Default Constructor
    public Vector3D()
    {
	x=0;
	y=0;
	z=0;
    }
    // Copy Constructor
    public Vector3D(Vector3D original)
    {
	this.setVector3D(original.getX(), original.getY(),original.getZ());
    }
    // Explicit Constructor
    public Vector3D(double x,double y,double z)
    {
	this.x=x;
	this.y=y;
	this.z=z;
    }

    //Methods
    public String toString()
    {
	return ("["+x+", "+y+", "+z+"]");
    }
    public double magnitude()
    {
	return Math.sqrt(x*x+y*y+z*z);
    }
    public double magnitudeSq()
    {
	return magnitude()*magnitude();
    }
    public Vector3D scale(double k)
    {
	return new Vector3D(getX()*k,getY()*k,getZ()*k);
    }
    public static Vector3D add(Vector3D a, Vector3D b)
    {
	return new Vector3D(a.getX()+b.getX(),a.getY()+b.getY(),a.getZ()+b.getZ());
    }
    public static Vector3D subtract(Vector3D a, Vector3D b)
    {
	return new Vector3D(a.getX()-b.getX(),a.getY()-b.getY(),a.getZ()-b.getZ());
    }   
    public static double dot(Vector3D a, Vector3D b)
    {
	return a.getX()*b.getX()+a.getY()*b.getY()+a.getZ()*b.getZ();
    }
    public static Vector3D cross(Vector3D a, Vector3D b)
    {
	return new Vector3D(a.getY()*b.getZ()-a.getZ()*b.getY(),
			    a.getZ()*b.getX()-a.getX()*b.getZ(),
			    a.getX()*b.getY()-a.getY()*b.getX());
    }
    public static boolean compare(Vector3D a, Vector3D b)
    {
	boolean test=false;
	if(a.getX()==b.getX() && a.getY()==b.getY() && a.getZ()==b.getZ())
	    {
		test=true;
	    }
	return test;
    }

    // Accessors and Mutators
    public void setVector3D(double xx, double yy, double zz)
    {
	x=xx;
	y=yy;
	z=zz;
    }
    public double getX()
    {
	return x;
    }
    public void setX(double xx)
    {
	x=xx;
    }
    public double getY()
    {
	return y;
    }
    public void setY(double yy)
    {
	y=yy;
    }
    public double getZ()
    {
	return z;
    }
    public void setZ(double zz)
    {
	z=zz;
    }
}
    