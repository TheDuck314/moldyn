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
    static final int imageWidth = 650;
    static final int imageHeight = 650;
    static double scale = 10;
    static final double width = imageWidth / scale;
    static final double height = imageHeight / scale;

    static Random rand = new Random(new Date().getTime());
    
    static final int N = 500;
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
    
    static int NcellX = 7;
    static int NcellY = 7;
    static double cellWidth = width / NcellX;
    static double cellHeight = height / NcellY;
    static int[][][] cellMembers = new int[NcellX][NcellY][1000];
    static int[][] cellPops = new int[NcellX][NcellY];
    static int[] cellNeighborOffsetsX = new int[] { 0, 1, 1, 0, -1 };
    static int[] cellNeighborOffsetsY = new int[] { 0, 0, 1, 1, 1 };

    static double pressure = 0.0001;
    static final double wallDensity = 1.0;
    static double[] wallCoords = new double[4];
    static double[] wallvs = new double[4];
    static double[] wallhalfvs = new double[4];
    static double[] wallas = new double[4];
    
    static boolean resampleVelocities = false;
    
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
                System.out.println("sup bitches");
                if (e.getKeyChar() == 'h') {
                    for (int i = 0; i < Main.N; ++i) {
                        Main.vxs[i] *= 1.1;
                        Main.vys[i] *= 1.1;
                    }
                    System.out.println("hotter");
                } else if (e.getKeyChar() == 'c') {
                    for (int i = 0; i < Main.N; ++i) {
                        Main.vxs[i] *= 0.9;
                        Main.vys[i] *= 0.9;
                    }               
                    System.out.println("colder");
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
        if (Math.sqrt(r2cutoff) >= cellWidth || Math.sqrt(r2cutoff) >= cellHeight) {
            System.out.println("cells too small");
            System.exit(-1);
        }
        
        System.out.format("NcellX = %d, NcellY = %d, cellWidth = %.2f, cellHeight = %.2f, rcutoff = %.2f\n", 
                NcellX, NcellY, cellWidth, cellHeight, Math.sqrt(r2cutoff));
        
        for (int i = 0; i < N; ++i) {
            //vxs[i] = 1*(-0.5 + Math.random());
            //vys[i] = 1*(-0.5*Math.random());

            boolean collision;
            do {
                xs[i] = /*5 + Math.random()*(width-10);*/ Math.random() * width;
                ys[i] = /*5 + Math.random()*(height-10);*/ Math.random() * height;
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

        wallCoords[0] = 0;
        wallCoords[1] = width;
        wallCoords[2] = 0;
        wallCoords[3] = height;

        double dt = 0.005;
        int iter = 0;
        
        long startTime = System.currentTimeMillis();
        
        while(true) {            
            if (iter % 1 == 0) drawAtoms();
           // System.in.read(); System.in.read();
            
            // compute half-step velocities
            for (int i = 0; i < N; ++i) {
                halfvxs[i] = vxs[i] + 0.5 * dt * axs[i];
                halfvys[i] = vys[i] + 0.5 * dt * ays[i];
                //System.out.format("halfvs[%d] = (%f,%f)\n", i, halfvxs[i], halfvys[i]);
            }
            for (int i = 0; i < 4; ++i) {
                wallhalfvs[i] = wallvs[i] + 0.5 * dt * wallas[i];
            }

            // compute new positions
            for (int i = 0; i < N; ++i) {
                xs[i] += dt * halfvxs[i];
                ys[i] += dt * halfvys[i];
//                if (xs[i] < 0) xs[i] += width;
//                if (ys[i] < 0) ys[i] += height;
//                if (xs[i] > width) xs[i] -= width;
//                if (ys[i] > height) ys[i] -= height;
                //System.out.format("x[%d] = (%f,%f)\n", i, xs[i], ys[i]);
            }
            for (int i = 0; i < 4; ++i) {
                wallCoords[i] += dt * wallhalfvs[i];
            }
            
            // compute new accelerations
            for (int i = 0; i < N; ++i) {
                axs[i] = 0;
                ays[i] = 0;
                checkaxs[i] = 0;
                checkays[i] = 0;
            }            
            wallas[0] = pressure / wallDensity;
            wallas[1] = -pressure / wallDensity;
            wallas[2] = pressure / wallDensity;
            wallas[3] = -pressure / wallDensity;
            //int zeroforcecountcheck = 0;
            /*for (int i = 0; i < N; ++i) {
                for (int j = 0; j < i; ++j) {
                    double dx = xs[i] - xs[j];
                    double dy = ys[i] - ys[j];
//                    if (dx > width/2) dx -= width;
//                    if (dy > height/2) dy -= height;
//                    if (dx < -width/2) dx += width;
//                    if (dy < -height/2) dy += height;
                    double r2 = dx*dx + dy*dy;
                    if (r2 <= r2cutoff) {
                        double rm6 = 1/(r2*r2*r2);
                        double Foverr = 12 * rm6 * (rm6 - 1);
                        double Fx = Foverr * dx;
                        double Fy = Foverr * dy;
                        checkaxs[i] += Fx;
                        checkays[i] += Fy; 
                        checkaxs[j] -= Fx;
                        checkays[j] -= Fy;
//                        if (i == 0) {
//                            System.out.format("check: %d <-> %d        ;    (%.4f, %.4f) <-> (%.4f, %.4f)    ; F(0) = (%.4f, %.4f)\n", i, j, xs[i], ys[i], xs[j], ys[j], Fx, Fy);
//                            zeroforcecountcheck += 1;
//                        }
                    }
                }
                //System.out.format("as[%d] = (%f,%f)\n", i, axs[i], ays[i]);
            }*/
            /*for (int i = 0; i < N; ++i) {
                axs[i] = checkaxs[i];
                ays[i] = checkays[i];
            }*/
            // put atoms in cells
            for (int cellX = 0; cellX < NcellX; ++cellX) {
                for(int cellY = 0; cellY < NcellY; ++cellY) {
                    cellPops[cellX][cellY] = 0;
                }
            }
            for (int i = 0; i < N; ++i) {
                int cellX = (int)(xs[i] / cellWidth);
                int cellY = (int)(ys[i] / cellHeight);
                if (cellX < 0) cellX = 0;
                if (cellY < 0) cellY = 0;
                if (cellX >= NcellX) cellX = NcellX - 1;
                if (cellY >= NcellY) cellY = NcellY - 1;
                cellMembers[cellX][cellY][cellPops[cellX][cellY]++] = i;
            }
            // compute forces between atoms
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
                                    double rm6 = 1/(r2*r2*r2);
                                    double Foverr = 12 * rm6 * (rm6 - 1);
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
            double wallStiffness = 5;
            //double buffer = 2;
            for (int i = 0; i < N; ++i) {
                //if (xs[i] < buffer) axs[i] += wallStiffness * (buffer - xs[i]);
                //if (ys[i] < buffer) ays[i] += wallStiffness * (buffer - ys[i]);
                //if (xs[i] > width-buffer) axs[i] -= wallStiffness * (xs[i] - width + buffer);
                //if (ys[i] > height-buffer) ays[i] -= wallStiffness * (ys[i] - height + buffer);
                if (xs[i] < wallCoords[0]) {
                    double F = wallStiffness * (wallCoords[0] - xs[i]);
                    axs[i] += F;
                    wallas[0] -= F / (wallDensity * height);
                }
                if (xs[i] > wallCoords[1]) {
                    double F = wallStiffness * (xs[i] - wallCoords[1]);
                    axs[i] -= F;
                    wallas[1] += F / (wallDensity * height);
                }
                if (ys[i] < wallCoords[2]) {
                    double F = wallStiffness * (wallCoords[2] - ys[i]);
                    ays[i] += F;
                    wallas[2] -= F / (wallDensity * width);
                }
                if (ys[i] > wallCoords[3]) {
                    double F = wallStiffness * (ys[i] - wallCoords[3]);
                    ays[i] -= F;
                    wallas[3] += F / (wallDensity * width);
                }
            }
//            System.out.println("zeroforcecount = " + zeroforcecount);
            /*for (int i = 0; i < N; ++i) {
                if (Math.abs(axs[i] - checkaxs[i]) > 1e-8) {
                    System.out.println("error in axs[" + i + "]; ax = " + axs[i] + ", checkax = " + checkaxs[i]);
                    System.exit(-1);
                }
                if (Math.abs(ays[i] - checkays[i]) > 1e-8) {
                    System.out.println("error in ays[" + i + "]");
                    System.exit(-1);
                }
            }*/
            
            // compute new velocities
            for (int i = 0; i < N; ++i) {
                vxs[i] = halfvxs[i] + 0.5 * dt * axs[i];
                vys[i] = halfvys[i] + 0.5 * dt * ays[i];
                //System.out.format("vs[%d] = (%f,%f)\n", i, vxs[i], vys[i]);
            }
            for (int i = 0; i < 4; ++i) {
                wallvs[i] = wallhalfvs[i] + 0.5 * dt * wallas[i];
            }
            if (resampleVelocities) {
                sampleVelocitiesFromMaxwellBoltzmann(1.0);
                resampleVelocities = false;
            }
            
            double comX = 0;
            double comY = 0;
            for (int i = 0; i < N; ++i) {
                comX += xs[i];
                comY += ys[i];
            }
            comX /= N;
            comY /= N;
            double shiftX = width*(10/scale)/2 - comX;
            double shiftY = height*(10/scale)/2 - comY;
            for (int i = 0; i < N; ++i) {
                xs[i] += shiftX;
                ys[i] += shiftY;
            }
            wallCoords[0] += shiftX;
            wallCoords[1] += shiftX;
            wallCoords[2] += shiftY;
            wallCoords[3] += shiftY;
            
            iter += 1;
            //if (iter % 10000 == 0) {
            //    sampleVelocitiesFromMaxwellBoltzmann(0.8);
            //}
            
            if (iter % 100 == 0) System.out.format("steps/sec = %.2f\n", (1000.0 * iter) / (System.currentTimeMillis() - startTime));
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
            
            g.fillOval((int)(scale * xs[i]-0.5*diameter), (int)(scale * ys[i]-0.5*diameter), diameter, diameter);
        }
        
        g.setColor(Color.white);
        g.drawLine((int)(scale * wallCoords[0]), 0, (int)(scale * wallCoords[0]), imageHeight);
        g.drawLine((int)(scale * wallCoords[1]), 0, (int)(scale * wallCoords[1]), imageHeight);
        g.drawLine(0, (int)(scale * wallCoords[2]), imageWidth, (int)(scale * wallCoords[2]));
        g.drawLine(0, (int)(scale * wallCoords[3]), imageWidth, (int)(scale * wallCoords[3]));
        
        g.setColor(Color.white);
        g.drawString(String.format("temp = %.4f", computeTemperature()), 10, 10);
        g.drawString(String.format("pressure = %.8f", pressure), 10, 30);
        
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
