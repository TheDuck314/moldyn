import java.awt.Color;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.WindowConstants;


public class Main {
    static JFrame frame;
    static BufferedImage image;

    static double scale = 10;
    static int N = 5000;
    static double T = 0.45;
    static double P = 0.03;

    public static void main(String[] args) throws IOException {
        
        final Sim sim = new MonatomicLennardJonesSim(N, T, P);
        
        frame = new JFrame("moldyn");
        frame.setSize(600, 600);
        frame.setVisible(true);
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        
        frame.addKeyListener(new KeyListener () {
            public void keyPressed(KeyEvent e) {}
            public void keyReleased(KeyEvent e) {}
            public void keyTyped(KeyEvent e) {
                if (e.getKeyChar() == 'h') {
                    T *= 1.01;
                    sim.SetT(T);
                } else if (e.getKeyChar() == 'c') {
                    T /= 1.01;
                    sim.SetT(T);
                } else if (e.getKeyChar() == 'p') {
                    P *= 1.03;
                    sim.SetP(P);
                } else if (e.getKeyChar() == 'o') {
                    P /= 1.03;
                    sim.SetP(P);
                } else if (e.getKeyChar() == '=') {
                    scale *= 2;
                } else if (e.getKeyChar() == '-') {
                    scale *= 0.5;
                } else if (e.getKeyChar() == 't') {
                    sim.ResampleVelocities();
                }
            }
        });
        
        double dt = 0.02; //0.005;
        
        double targetFrameRate = 30;
        double millisPerStep = 1;
        int stepsSinceFrame = Integer.MAX_VALUE;
        while(true) {                 
            long stepStartTime = System.currentTimeMillis();

            double numSkip = 1000.0 / (millisPerStep * targetFrameRate);
            if (stepsSinceFrame >= numSkip) {
                Draw(sim);
                stepsSinceFrame = 0;
            } else {
                stepsSinceFrame += 1;
            }
            
            //System.in.read(); System.in.read();

                      
            sim.ComputeHalfVelocities(dt);
            sim.ComputeNewPositions(dt);
            sim.ComputeAcclerations(dt);
            sim.ComputeNewVelocities(dt);
            sim.EndStep(dt);
            
            long stepMillis = System.currentTimeMillis() - stepStartTime;
            millisPerStep = 0.9999 * millisPerStep + 0.0001 * stepMillis;
        }
    }
    
    static void Draw(Sim sim)
    {
        int imageWidth = frame.getWidth() - 20;
        int imageHeight = frame.getHeight() - 50;
        if (image == null || image.getWidth() != imageWidth || image.getHeight() != imageHeight) {
            image = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_INT_RGB);
        }
        
        sim.Draw(image.getGraphics(), imageWidth, imageHeight, scale);
        
        frame.getGraphics().drawImage(image, 9, 32, frame);
    }
    

    

}
