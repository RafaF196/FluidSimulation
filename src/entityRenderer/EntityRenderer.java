package entityRenderer;

import java.util.List;

import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.GL13;
import org.lwjgl.opengl.GL20;
import org.lwjgl.opengl.GL30;
import org.lwjgl.opengl.GL31;
import org.lwjgl.util.vector.Matrix4f;

import entities.Camera;
import entities.Entity;
import models.RawModel;
import models.TexturedModel;
import toolbox.MatrixOps;

public class EntityRenderer {

	private EntityShader entityShader;

	public EntityRenderer(Matrix4f projectionMatrix) {
		entityShader = new EntityShader();
		entityShader.start();
		entityShader.loadProjectionMatrix(projectionMatrix);
		entityShader.connectTextureUnits();
		entityShader.stop();
	}

	public void render(List<Entity> entities, Camera camera, Integer N) {
		entityShader.start();
		entityShader.loadViewMatrix(camera);
		entityShader.loadOffsetPosition(N);
		for (Entity entity : entities) {
			TexturedModel model = entity.getModel();
			bindModelVao(model);
			bindTexture(model);
			loadModelMatrix(entity);
			// used "glDrawElementsInstanced" for creating several copies of the same object
			// the number of copies will be defined by N, which is the number of elements in one row.
			// Hence, the total amount of copies will be N*N
			GL31.glDrawElementsInstanced(GL11.GL_TRIANGLES, model.getRawModel().getVertexCount(),
					GL11.GL_UNSIGNED_INT, 0, N*N);
			unbindVao();
		}
		entityShader.stop();
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

	private void loadModelMatrix(Entity entity) {
		Matrix4f transformationMatrix = MatrixOps.createTransformationMatrix(entity.getPosition(), 0, entity.getRotY(), 0,
				entity.getScale());
		entityShader.loadTransformationMatrix(transformationMatrix);
	}
	
	public void cleanUp() {
		entityShader.cleanUp();
	}

}
