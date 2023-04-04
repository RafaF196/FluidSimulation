package entities;

import org.lwjgl.input.Keyboard;
import org.lwjgl.input.Mouse;
import org.lwjgl.util.vector.Vector3f;

public class Camera {
	
	private Vector3f position;
	private float pitch;
	private float yaw;
	
	public Camera(Vector3f position, Float pitch, Float yaw) {
		this.position = position;
		this.pitch = pitch;
		this.yaw = yaw;
	}
	
	public void move(Float delta) {
		moveCamera(delta);
		calculatePitch();
		calculateYaw();
	}
	
	private void moveCamera(Float delta) {
		Vector3f actPos = this.position;
		Float step = 5 * delta;
		if (Keyboard.isKeyDown(Keyboard.KEY_D)) {
			this.position = new Vector3f(actPos.x + step, actPos.y, actPos.z);
		}
		if (Keyboard.isKeyDown(Keyboard.KEY_A)) {
			this.position = new Vector3f(actPos.x - step, actPos.y, actPos.z);
		}
		if (Keyboard.isKeyDown(Keyboard.KEY_W)) {
			this.position = new Vector3f(actPos.x, actPos.y, actPos.z - step);
		}
		if (Keyboard.isKeyDown(Keyboard.KEY_S)) {
			this.position = new Vector3f(actPos.x, actPos.y, actPos.z + step);
		}
		if (Keyboard.isKeyDown(Keyboard.KEY_Q)) {
			this.position = new Vector3f(actPos.x, actPos.y - step, actPos.z);
		}
		if (Keyboard.isKeyDown(Keyboard.KEY_E)) {
			this.position = new Vector3f(actPos.x, actPos.y + step, actPos.z);
		}
	}
	
	private void calculateYaw() {
		if (Mouse.isButtonDown(0)) {
			float angleChange = Mouse.getDX() * 0.2f;
			yaw += angleChange;
			yaw %= 360;
		}
	}
	
	private void calculatePitch() {
		if (Mouse.isButtonDown(0)) {
			float pitchChange = Mouse.getDY() * 0.2f;
			pitch -= pitchChange;
			if (pitch < -90){
				pitch = -90;
			} else if (pitch > 90){
				pitch = 90;
			}
		}
	}

	public void invertPitch() {
		this.pitch = -pitch;
	}

	public Vector3f getPosition() {
		return position;
	}

	public float getPitch() {
		return pitch;
	}

	public float getYaw() {
		return yaw;
	}

}
