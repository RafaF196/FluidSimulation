package collisions;

import org.lwjgl.util.vector.Vector3f;

public class SphereCol {
	
	private Vector3f center;
	private Float radius;
	
	public SphereCol(Vector3f center, Float radius) {
		this.center = center;
		this.radius = radius;
	}
	
	public Vector3f getCenter() {
		return center;
	}
	
	public Float getRadius() {
		return radius;
	}

}
