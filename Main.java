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


public class Main {
    static JFrame frame;
    static BufferedImage image;
    static final int imageWidth = 1000;
    static final int imageHeight = 1000;
    static double scale = 8;
    //static final double width = imageWidth / scale;
    //static final double height = imageHeight / scale;

    static Random rand = new Random(new Date().getTime());
    
    static final int N = 5000;
    static double[] xs = new double[N];
    static double[] ys = new double[N];
    static double[] vxs = new double[N];
    static double[] vys = new double[N];
    static double[] axs = new double[N];
    static double[] ays = new double[N];
    static double[] halfvxs = new double[N];
    static double[] halfvys = new double[N];
    static double[] checkaxs = new double[N];
    static double[] checkays = new double[N];
    
    static int NcellX;
    static int NcellY;
    static double cellWidth;
    static double cellHeight;
    static int[][][] cellMembers;
    static int[][] cellPops;
    static int[] cellNeighborOffsetsX = new int[] { 0, 1, 1, 0, -1 };
    static int[] cellNeighborOffsetsY = new int[] { 0, 0, 1, 1, 1 };

    static double pressure = 0.1;
    static double wallRadius;
    static double wallMass = 50;
    static double wallStiffness = 5;
    static double wallV;
    static double wallHalfV;
    static double wallA;
    
    static double targetT = 1.0;
    static double avgT;
    
    static double totalInteratomPE;
    static double totalAtomKE;
    static double totalWallPE;
    static double wallKE;
    static double pressureE;
    static double totalE;
    
    static boolean resampleVelocities = false;
    
    static double millisPerStep = 1.0;
    
    public static void main(String[] args) throws IOException {
        System.out.println("hello world!");
        
        image = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_INT_RGB);
        
        frame = new JFrame("moldyn");
        frame.setSize(imageWidth + 50, imageHeight + 50);
        frame.setVisible(true);
        
        frame.addKeyListener(new KeyListener () {
            public void keyPressed(KeyEvent e) {}
            public void keyReleased(KeyEvent e) {}
            public void keyTyped(KeyEvent e) {
                if (e.getKeyChar() == 'h') {
                    Main.targetT *= 1.01;
                } else if (e.getKeyChar() == 'c') {
                    Main.targetT *= 0.99;
                } else if (e.getKeyChar() == 'p') {
                    Main.pressure *= 1.1;
                } else if (e.getKeyChar() == 'o') {
                    Main.pressure *= 0.9;
                } else if (e.getKeyChar() == '=') {
                    Main.scale *= 2;
                } else if (e.getKeyChar() == '-') {
                    Main.scale *= 0.5;
                } else if (e.getKeyChar() == 't') {
                    Main.resampleVelocities = true;
                }
            }
        });
        
        /*for (int x = 0; x < 100; ++x) {
            for(int y = 0; y < 100; ++y) {
                image.setRGB(x, y, 0x12345678);
            }
        }
        Graphics g = image.getGraphics();
        g.setColor(Color.red);
        g.fillRect(50, 50, 100, 100);
        frame.getGraphics().drawImage(image, 0, 0, frame);*/
         
       
        double r2cutoff = 9;
        double rcutoff = Math.sqrt(r2cutoff);
        
        System.out.format("NcellX = %d, NcellY = %d, cellWidth = %.2f, cellHeight = %.2f, rcutoff = %.2f\n", 
                NcellX, NcellY, cellWidth, cellHeight, Math.sqrt(r2cutoff));
        
        wallRadius = Math.max(imageWidth, imageHeight) / (2*scale);
        
        for (int i = 0; i < N; ++i) {
            //vxs[i] = 1*(-0.5 + Math.random());
            //vys[i] = 1*(-0.5*Math.random());

            boolean collision;
            do {
                xs[i] = 2 * (-0.5 + Math.random()) * wallRadius;
                ys[i] = 2 * (-0.5 + Math.random()) * wallRadius;
                collision = false;
                for (int j = 0; j < i; ++j) {
                    double dx = xs[i] - xs[j];
                    double dy = ys[i] - ys[j];
//                    if (dx > width/2) dx -= width;
//                    if (dy > height/2) dy -= height;
//                    if (dx < -width/2) dx += width;
//                    if (dy < -height/2) dy += height;
                    if (Math.sqrt(dx*dx + dy*dy) < 1) {
                        collision = true;
                        break;
                    }
                }            
            } while(collision);
        }
        sampleVelocitiesFromMaxwellBoltzmann(0.8);
        wallV = 0;
        
/*        xs[0] = -1;
        xs[1] = 1;
        ys[0] = ys[1] = 0;
        vxs[0] = vys[0] = vxs[1] = vys[1] = 0;*/
        
        double dt = 0.005;
        int iter = 0;
        
        
        double targetFrameRate = 60;
        millisPerStep = 1;
        int stepsSinceFrame = Integer.MAX_VALUE;
        avgT = computeTemperature();
        while(true) {     
            long stepStartTime = System.currentTimeMillis();

            double numSkip = 1000.0 / (millisPerStep * targetFrameRate);
            if (stepsSinceFrame >= numSkip) {
                drawAtoms();
                stepsSinceFrame = 0;
            } else {
                stepsSinceFrame += 1;
            }
            
            // compute half-step velocities
            for (int i = 0; i < N; ++i) {
                halfvxs[i] = vxs[i] + 0.5 * dt * axs[i];
                halfvys[i] = vys[i] + 0.5 * dt * ays[i];
            }
            wallHalfV = wallV + 0.5 * dt * wallA;

            // compute new positions
            for (int i = 0; i < N; ++i) {
                xs[i] += dt * halfvxs[i];
                ys[i] += dt * halfvys[i];
            }
            wallRadius += dt * wallHalfV;
            
            // compute new accelerations
            for (int i = 0; i < N; ++i) {
                axs[i] = 0;
                ays[i] = 0;
                checkaxs[i] = 0;
                checkays[i] = 0;
            }            
            double wallArea = 8*wallRadius;
            double wallForce = -pressure * wallArea;
            wallA = wallForce / wallMass;
            // the force law for the wall radius is d^2 R/dt^2 = -8*R*P/M
            // This can be rewritten dS/dt = -8*R*P where S = M*Rdot is the wall momentum
            // The corresponding Hamiltonian is H = S^2/(2M) + 4*R^2*P, that is H = S^2/2M + PV
            // where V = (2R)^2 is the volume.
            double volume = 4*wallRadius*wallRadius;
            pressureE = pressure * volume;

            // decide on cell size
            if (iter % 10 == 0) {
                NcellX = (int)Math.max(1, Math.min(15, (2*wallRadius / rcutoff)));
                NcellY = (int)Math.max(1, Math.min(15, (2*wallRadius / rcutoff)));
                cellWidth = 2 * wallRadius / NcellX;
                cellHeight = 2 * wallRadius / NcellY;
                if (cellWidth <= rcutoff || cellHeight <= rcutoff) {
                    System.out.println("cells too small");
                    System.exit(-1);;
                } else {
                    //System.out.format("NcellX = %d, NcellY = %d, cells are %.2f x %.2f\n", NcellX, NcellY, cellWidth, cellHeight);
                }
                cellMembers = new int[NcellX][NcellY][1000];
                cellPops = new int[NcellX][NcellY];
            }            

            
            // put atoms in cells
            for (int cellX = 0; cellX < NcellX; ++cellX) {
                for(int cellY = 0; cellY < NcellY; ++cellY) {
                    cellPops[cellX][cellY] = 0;
                }
            }
            for (int i = 0; i < N; ++i) {
                int cellX = (int)((wallRadius + xs[i]) / cellWidth);
                int cellY = (int)((wallRadius + ys[i]) / cellHeight);
                if (cellX < 0) cellX = 0;
                if (cellY < 0) cellY = 0;
                if (cellX >= NcellX) cellX = NcellX - 1;
                if (cellY >= NcellY) cellY = NcellY - 1;
                cellMembers[cellX][cellY][cellPops[cellX][cellY]++] = i;
            }
            // compute forces between atoms
            totalInteratomPE = 0;
            for (int cellX = 0; cellX < NcellX; ++cellX) {
                for (int cellY = 0; cellY < NcellY; ++cellY) {
                    int[] thisCellMembers = cellMembers[cellX][cellY];
                    int thisCellPop = cellPops[cellX][cellY];
                    for (int off = 0; off < 5; ++off) {
                        int neighborCellX = (cellX + cellNeighborOffsetsX[off] + NcellX) % NcellX;
                        int neighborCellY = (cellY + cellNeighborOffsetsY[off] + NcellY) % NcellY;
                        int[] neighborCellMembers = cellMembers[neighborCellX][neighborCellY];
                        int neighborCellPop = cellPops[neighborCellX][neighborCellY];
                        for (int thisAtom = 0; thisAtom < thisCellPop; ++thisAtom) {
                            for (int neighborAtom = 0; neighborAtom < neighborCellPop; ++neighborAtom) {
                                int i = thisCellMembers[thisAtom];
                                int j = neighborCellMembers[neighborAtom];
                                if (off == 0 && neighborAtom >= thisAtom) break; // avoid double-counting forces within a cell
                                
                                double dx = xs[i] - xs[j];
                                double dy = ys[i] - ys[j];
//                                if (dx > width/2) dx -= width;
//                                if (dy > height/2) dy -= height;
//                                if (dx < -width/2) dx += width;
//                                if (dy < -height/2) dy += height;
                                double r2 = dx*dx + dy*dy;
                                if (r2 <= r2cutoff) {
//                                    if (i == 0 || j == 0) {
//                                        System.out.format("cells: %d <-> %d        ;    (%.4f, %.4f) <-> (%.4f, %.4f)\n", i, j, xs[i], ys[i], xs[j], ys[j]);
//                                        zeroforcecount += 1;
//                                    }
                                    // V = 1/r^12 - 2/r^6
                                    // F = -dV/dr = 12 (1/r^13 - 1/r^7)
                                    // F/r = 12 * 1/r^2 * (1/r^12 - 1/r^6) = 12 (1/r^2) (1/r^6) (1/r^6 - 1)
                                    double rm6 = 1/(r2*r2*r2);
                                    double Foverr = 12 * rm6 * (rm6 - 1) / r2;
                                    totalInteratomPE += rm6 * (rm6 - 2);
                                    double Fx = Foverr * dx;
                                    double Fy = Foverr * dy;
                                    axs[i] += Fx;
                                    ays[i] += Fy; 
                                    axs[j] -= Fx;
                                    ays[j] -= Fy;
                                }
                            }
                        }
                    }                    
                }
            }
            // wall forces
            totalWallPE = 0;
            for (int i = 0; i < N; ++i) {
                if (xs[i] < -wallRadius) {
                    double d = (xs[i] + wallRadius);
                    double F = -wallStiffness * d;
                    axs[i] += F;
                    wallA += F / wallMass;
                    totalWallPE += 0.5 * wallStiffness * d * d;
                }
                if (xs[i] > wallRadius) {
                    double d = (xs[i] - wallRadius);
                    double F = wallStiffness * d;
                    axs[i] -= F;
                    wallA += F / wallMass;
                    totalWallPE += 0.5 * wallStiffness * d * d;
                }
                if (ys[i] < -wallRadius) {
                    double d = (ys[i] + wallRadius);
                    double F = -wallStiffness * d;
                    ays[i] += F;
                    wallA += F / wallMass;
                    totalWallPE += 0.5 * wallStiffness * d * d;
                }
                if (ys[i] > wallRadius) {
                    double d = (ys[i] - wallRadius);
                    double F = wallStiffness * d;
                    ays[i] -= F;
                    wallA += F / wallMass;
                    totalWallPE += 0.5 * wallStiffness * d * d;
                }
            }
            
            // compute new velocities
            totalAtomKE = 0;
            for (int i = 0; i < N; ++i) {
                vxs[i] = halfvxs[i] + 0.5 * dt * axs[i];
                vys[i] = halfvys[i] + 0.5 * dt * ays[i];
                totalAtomKE += 0.5 * (vxs[i] * vxs[i] + vys[i] * vys[i]);
            }
            wallV = wallHalfV + 0.5 * dt * wallA;
            if (resampleVelocities || iter % 2000 == 0) {
                sampleVelocitiesFromMaxwellBoltzmann(targetT);
                resampleVelocities = false;
            }
            wallKE = 0.5 * wallMass * wallV * wallV;
            
            totalE = totalAtomKE + totalInteratomPE + totalWallPE + wallKE + pressureE;
            
            iter += 1;
            
            avgT = 0.999 * avgT + 0.001 * computeTemperature();
            
            long stepMillis = System.currentTimeMillis() - stepStartTime;
            millisPerStep = 0.9999 * millisPerStep + 0.0001 * stepMillis;
        }
    }
    
    static void drawAtoms()
    {
        int diameter = (int)(scale);
        if (diameter < 2) diameter = 2;
        double refspeed = 2;
        
        Graphics g = image.getGraphics();
        g.setColor(Color.black);
        g.fillRect(0, 0, imageWidth, imageHeight);
        
        g.setColor(Color.blue);
        for (int i = 0; i < xs.length; ++i) {            
            double speed = Math.sqrt(vxs[i] * vxs[i] + vys[i] * vys[i]);
            double red = 0, green = 0, blue = 0;
            if (speed < refspeed) {
                blue = 1 - speed/refspeed;
                red = speed/refspeed;
            } else if (speed < 2*refspeed) {
                red = 1;
                green = speed/refspeed - 1;
            } else if (speed < 3*refspeed) {
                green = 1;
                red = speed/refspeed - 2;
            } else {
                red = 1;
                green = 1;
            }
            g.setColor(new Color((float)red, (float)green, (float)blue));
            
            g.fillOval(imageWidth/2 + (int)(scale * xs[i]-0.5*diameter), imageHeight/2 + (int)(scale * ys[i]-0.5*diameter), diameter, diameter);
        }
        
        g.setColor(Color.white);
        g.drawRect(imageWidth/2 - (int)(scale * wallRadius), imageHeight/2 - (int)(scale * wallRadius), (int)(2*scale*wallRadius), (int)(2*scale*wallRadius));
        
        g.setColor(Color.white);
        g.drawString(String.format("targetT = %.4f, temp = %.4f", targetT, avgT), 10, 10);
        g.drawString(String.format("pressure = %.8f", pressure), 10, 30);
        g.drawString(String.format("step = %.2f ms", millisPerStep), 10, 50);
        g.drawString(String.format("totalE = %.2f", totalE), 10, 70);
        g.drawString(String.format("totalAtomKE = %.2f", totalAtomKE), 10, 90);
        g.drawString(String.format("totalInteratomPE = %.2f", totalInteratomPE), 10, 110);
        g.drawString(String.format("totalWallPE = %.2f", totalWallPE), 10, 130);
        g.drawString(String.format("wallKE = %.2f", wallKE), 10, 150);
        g.drawString(String.format("pressureE = %.2f", pressureE), 10, 170);
        
        frame.getGraphics().drawImage(image, 9, 32, frame);
    }
    
    static double computeTemperature() 
    {
        double meanVx2 = 0;
        double meanVy2 = 0;
        for (int i = 0; i < vxs.length; ++i) {
            meanVx2 += vxs[i] * vxs[i];
            meanVy2 += vys[i] * vys[i];
        }
        meanVx2 /= vxs.length;
        meanVy2 /= vxs.length;
        return (meanVx2 + meanVy2)/2;
    }
    
    // Maxwell-Boltzmann distribution for each velocity component is gaussian with mean 0
    // and 
    static void sampleVelocitiesFromMaxwellBoltzmann(double T)
    {
        double std = Math.sqrt(T);
        for (int i = 0; i < vxs.length; ++i) {
            vxs[i] = std * rand.nextGaussian();
            vys[i] = std * rand.nextGaussian();
        }
        System.out.format("set velocities for T = %f; computed temp is %.4f\n", T, computeTemperature());
    }
}
