package render;

import java.util.List;

import org.lwjgl.opengl.Display;
import org.lwjgl.opengl.GL11;
import org.lwjgl.util.vector.Matrix4f;

import entities.Camera;
import entities.Entity;
import entityRenderer.EntityRenderer;
import particles.ParticleRenderer;
import particles.ParticleSystem;

public class MasterRenderer {
	
	private static final float FOV = 70;
	private static final float NEAR_PLANE = 0.1f;
	private static final float FAR_PLANE = 1000;
	
	private Matrix4f projectionMatrix;
	
	private EntityRenderer entityRenderer;
	private ParticleRenderer particleRenderer;
	
	public MasterRenderer(Loader loader) {
		this.projectionMatrix = createProjectionMatrix();
		this.entityRenderer = new EntityRenderer(projectionMatrix);
		this.particleRenderer = new ParticleRenderer(projectionMatrix);
	}
	
	public void renderScene(List<Entity> entities, ParticleSystem ps, Camera camera, Integer N) {
		prepare();
		entityRenderer.render(entities, camera, N);
		particleRenderer.render(ps, camera);
	}
	
	public void cleanUp() {
		entityRenderer.cleanUp();
		particleRenderer.cleanUp();
	}
	
	private static Matrix4f createProjectionMatrix() {
		Matrix4f projectionMatrix = new Matrix4f();
		float aspectRatio = (float) Display.getWidth() / (float) Display.getHeight();
		float y_scale = (float) ((1f / Math.tan(Math.toRadians(FOV / 2f))));
		float x_scale = y_scale / aspectRatio;
		float frustum_length = FAR_PLANE - NEAR_PLANE;

		projectionMatrix.m00 = x_scale;
		projectionMatrix.m11 = y_scale;
		projectionMatrix.m22 = -((FAR_PLANE + NEAR_PLANE) / frustum_length);
		projectionMatrix.m23 = -1;
		projectionMatrix.m32 = -((2 * NEAR_PLANE * FAR_PLANE) / frustum_length);
		projectionMatrix.m33 = 0;
		return projectionMatrix;
	}
	
	private static void prepare() {
		GL11.glEnable(GL11.GL_DEPTH_TEST);
		GL11.glClear(GL11.GL_COLOR_BUFFER_BIT|GL11.GL_DEPTH_BUFFER_BIT);
		GL11.glEnable(GL11.GL_CULL_FACE);
		GL11.glCullFace(GL11.GL_BACK);
	}

}
