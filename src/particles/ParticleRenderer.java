package particles;

import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.GL13;
import org.lwjgl.opengl.GL20;
import org.lwjgl.opengl.GL30;
import org.lwjgl.opengl.GL31;
import org.lwjgl.util.vector.Matrix4f;

import entities.Camera;
import models.RawModel;
import models.TexturedModel;
import toolbox.MatrixOps;

public class ParticleRenderer {

	private ParticleShader particleShader;

	public ParticleRenderer(Matrix4f projectionMatrix) {
		particleShader = new ParticleShader();
		particleShader.start();
		particleShader.loadProjectionMatrix(projectionMatrix);
		particleShader.connectTextureUnits();
		particleShader.stop();
	}

	public void render(ParticleSystem ps, Camera camera) {
		particleShader.start();
		particleShader.loadViewMatrix(camera);
		particleShader.loadOffsetPosition(ps.getOffsets());
		
		TexturedModel model = ps.getModel();
		bindModelVao(model);
		bindTexture(model);
		loadModelMatrix(ps);
		GL31.glDrawElementsInstanced(GL11.GL_TRIANGLES, model.getRawModel().getVertexCount(),
				GL11.GL_UNSIGNED_INT, 0, ps.particles.size());
		unbindVao();
			
		particleShader.stop();
	}

	// binds the Vertex Array Object of a model and enables the vertex attributes stored in positions
	// 0 (vertex positions), 1 (texture coordinates) and 2 (normal vectors)
	private void bindModelVao(TexturedModel model) {
		RawModel rawModel = model.getRawModel();
		GL30.glBindVertexArray(rawModel.getVaoID());
		GL20.glEnableVertexAttribArray(0);
		GL20.glEnableVertexAttribArray(1);
		GL20.glEnableVertexAttribArray(2);
	}

	// unbinds the Vertex Array Object
	private void unbindVao() {
		GL20.glDisableVertexAttribArray(0);
		GL20.glDisableVertexAttribArray(1);
		GL20.glDisableVertexAttribArray(2);
		GL30.glBindVertexArray(0);
	}

	// binds the texture of the model
	private void bindTexture(TexturedModel model) {
		GL13.glActiveTexture(GL13.GL_TEXTURE0);
		GL11.glBindTexture(GL11.GL_TEXTURE_2D, model.getTexture().getID());
	}

	private void loadModelMatrix(ParticleSystem ps) {
		Matrix4f transformationMatrix = MatrixOps.createTransformationMatrix(ps.getPosition(), 0, 0, 0, 0.05f);
		particleShader.loadTransformationMatrix(transformationMatrix);
	}
	
	public void cleanUp() {
		particleShader.cleanUp();
	}

}
