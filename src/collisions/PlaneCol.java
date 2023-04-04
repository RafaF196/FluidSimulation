package collisions;

import org.lwjgl.util.vector.Vector3f;

public class PlaneCol {
	
	private Vector3f normal;
	private Float value;
	
	public PlaneCol(Vector3f orientation, Float value) {
		this.normal = orientation;
		this.value = value;
	}
	
	public Vector3f getNormal() {
		return normal;
	}
	
	public Float getValue() {
		return value;
	}

}
