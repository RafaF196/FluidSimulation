package render;

import org.lwjgl.LWJGLException;
import org.lwjgl.Sys;
import org.lwjgl.opengl.Display;
import org.lwjgl.opengl.DisplayMode;
import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.PixelFormat;

public class DisplayManager {
	
	private static final int WIDTH = 1280;
	private static final int HEIGHT = 640;
	private static final int FPS_CAP = 100;
	
	private static long lastFrameTime; // the time when the last frame was rendered
	private static float delta; // time elapsed between last frame and the current frame (in seconds)
	
	public static void createDisplay() {			
		try {
			Display.setDisplayMode(new DisplayMode(WIDTH,HEIGHT));
			Display.setInitialBackground(0.75f, 0.75f, 0.8f);
			Display.setTitle("Fluid Simulation");
			Display.create(new PixelFormat().withDepthBits(24).withSamples(4));
		} catch (LWJGLException e) {
			e.printStackTrace();
		}
		
		GL11.glViewport(0, 0, WIDTH, HEIGHT);
	}
	
	// updates the Display (called every frame)
	public static void updateDisplay() {
		Display.sync(FPS_CAP);
		Display.update();
		long currentFrameTime = getCurrentTime();
		delta = (currentFrameTime - lastFrameTime)/1000f;
		lastFrameTime = currentFrameTime;
	}
	
	// returns the current time in milliseconds
	private static long getCurrentTime(){
		return Sys.getTime()*1000/Sys.getTimerResolution();
	}
	
	public static float getFrameTimeSeconds(){
		return delta;
	}
	
	public static void closeDisplay() {
		Display.destroy();
	}

}
