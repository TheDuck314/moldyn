import java.awt.Color;
import java.awt.Graphics;
import java.util.Date;
import java.util.Random;


public class DipoleSim implements Sim {
    static Random rand = new Random(new Date().getTime());
    
    int N;
    double[] xs;
    double[] ys;
    double[] vxs;
    double[] vys;
    double[] halfvxs;
    double[] halfvys;
    double[] axs;
    double[] ays;    
    
    double[] thetas;
    double[] omegas;
    double[] halfomegas;
    double[] torques;
    double momentOfInertia = 2;

    static final double dipoleRadius = 0.04;
    static final double dipoleStrength = 2.0;
    
    int[][][] cellMembers;
    int[][] cellPops;
    int[] cellNeighborOffsetsX = new int[] { 0, 1, 1, 0, -1 };
    int[] cellNeighborOffsetsY = new int[] { 0, 0, 1, 1, 1 };

    double pressure;
    double wallMass = 100;
    double wallStiffness = 5;
    double wallRadius;
    double wallV;
    double wallHalfV;
    double wallA;
    
    double targetT = 0.10;
    double avgT;
    
    double totalInteratomPE;
    double totalAtomKE;
    double totalRotationalKE;
    double totalWallPE;
    double wallKE;
    double pressureE;
    double totalE;

    double r2cutoff = 9;
    double rcutoff = Math.sqrt(r2cutoff);
    
    int iter = 0;
    double simTime = 0;
    
    boolean resampleVelocities = false;
    
    public DipoleSim(int N, double T, double P) {
        this.N = N;
        
        xs = new double[N];
        ys = new double[N];
        vxs = new double[N];
        vys = new double[N];
        halfvxs = new double[N];
        halfvys = new double[N];
        axs = new double[N];
        ays = new double[N];
        
        thetas = new double[N];
        omegas = new double[N];
        halfomegas = new double[N];
        torques = new double[N];

        
        pressure = P;
        targetT = T;
        
        wallRadius = 10 + Math.sqrt(N);
        
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
                thetas[i] = 2 * Math.PI * Math.random();
                omegas[i] = 0;
            } while(collision);
        }
        
        /*
        int rowLen = (int)Math.ceil(Math.sqrt(Math.sqrt(3)*N/2));
        double xSpacing = 1.0;
        double ySpacing = xSpacing * Math.sqrt(3)/2;
        double cornerOffset = 0.5 * xSpacing * rowLen;
        double oddRowOffset = xSpacing / 2;
        int row = 0;
        int col = 0;
        for (int i = 0; i < N; ++i) {
            xs[i] = col * xSpacing - cornerOffset + (row % 2 == 0 ? 0 : oddRowOffset);
            ys[i] = row * ySpacing - cornerOffset;
            col += 1;
            if (col == rowLen) {
                col = 0;
                row += 1;
            }
        }        
        wallRadius = cornerOffset + xSpacing;
        */
        
        wallV = 0;
        sampleVelocitiesFromMaxwellBoltzmann(targetT);
        
        
        /*xs[0] = -1.0;
        ys[0] = 0.0;
        xs[1] = +1.0;
        ys[1] = 0.0;
        vxs[0] = 0;
        vys[0] = 0;
        vxs[1] = 0;
        vys[1] = 0;
        thetas[0] = 0;
        thetas[1] = Math.PI/4;
        omegas[0] = 0;
        omegas[1] = 0;*/
        
        avgT = computeTemperature();
        
        int maxNcell = (int)Math.sqrt(N);
        cellMembers = new int[maxNcell][maxNcell][1000];
        cellPops = new int[maxNcell][maxNcell];
    }

    public void SetT(double T) {
        targetT = T;
    }
    
    public void SetP(double P) {
        pressure = P;
    }
    
    public void ResampleVelocities() {
        resampleVelocities = true;
    }
    
    public void ComputeHalfVelocities(double dt) {
        for (int i = 0; i < N; ++i) {
            halfvxs[i] = vxs[i] + 0.5 * dt * axs[i];
            halfvys[i] = vys[i] + 0.5 * dt * ays[i];
            halfomegas[i] = omegas[i] + 0.5 * dt * torques[i] / momentOfInertia;
        }
        wallHalfV = wallV + 0.5 * dt * wallA;
    }

    public void ComputeNewPositions(double dt) {
        // compute new positions
        for (int i = 0; i < N; ++i) {
            xs[i] += dt * halfvxs[i];
            ys[i] += dt * halfvys[i];
            thetas[i] += dt * halfomegas[i];
            if (thetas[i] < 0) thetas[i] += 2 * Math.PI;
            if (thetas[i] > 2 * Math.PI) thetas[i] -= 2 * Math.PI;
        }
        wallRadius += dt * wallHalfV;
    }
    
    void DoChargePair(int i, int j, int iSign, int jSign) {        
        double iChargeX = xs[i] + iSign * dipoleRadius * Math.cos(thetas[i]);
        double iChargeY = ys[i] + iSign * dipoleRadius * Math.sin(thetas[i]);
        double jChargeX = xs[j] + jSign * dipoleRadius * Math.cos(thetas[j]);
        double jChargeY = ys[j] + jSign * dipoleRadius * Math.sin(thetas[j]);
        
        double dx = iChargeX - jChargeX;
        double dy = iChargeY - jChargeY;
        double r2 = dx * dx + dy * dy;
        double rm2 = 1 / r2;
        double rm6 = rm2 * rm2 * rm2;
        double FoverR = 12 * dipoleStrength * iSign * jSign * rm6 * rm2;
        totalInteratomPE += 2 * dipoleStrength * iSign * jSign * rm6;
        double Fx = FoverR * dx;
        double Fy = FoverR * dy;
        
        axs[i] += Fx;
        ays[i] += Fy;
        axs[j] -= Fx;
        ays[j] -= Fy;
        
        // Torque = F cross r
        torques[i] += -(Fx * (iChargeY - ys[i]) - Fy * (iChargeX - xs[i]));
        torques[j] +=  (Fx * (jChargeY - ys[j]) - Fy * (jChargeX - xs[j]));
        
        //System.out.format("%d<->%d  (%d,%d) force, r=%f, F=%f\n", i, j, iSign, jSign, Math.sqrt(r2), FoverR*Math.sqrt(r2));
    }
    
    public void ComputeAcclerations(double dt) {
        // compute new accelerations
        for (int i = 0; i < N; ++i) {
            axs[i] = 0;
            ays[i] = 0;
            torques[i] = 0;
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
        int Ncell = Math.min((int)(2*wallRadius / rcutoff)-1, (int)Math.sqrt(N));
        if (Ncell < 1) Ncell = 1;
        double cellSize = 2 * wallRadius / Ncell;
        if (Ncell > 1 && cellSize <= rcutoff) {
            System.out.println("cells too small");
            System.exit(-1);
        }
        
        // put atoms in cells
        for (int cellX = 0; cellX < Ncell; ++cellX) {
            for(int cellY = 0; cellY < Ncell; ++cellY) {
                cellPops[cellX][cellY] = 0;
            }
        }
        for (int i = 0; i < N; ++i) {
            int cellX = (int)((wallRadius + xs[i]) / cellSize);
            int cellY = (int)((wallRadius + ys[i]) / cellSize);
            if (cellX < 0) cellX = 0;
            if (cellY < 0) cellY = 0;
            if (cellX >= Ncell) cellX = Ncell - 1;
            if (cellY >= Ncell) cellY = Ncell - 1;
            cellMembers[cellX][cellY][cellPops[cellX][cellY]++] = i;
        }
        // compute forces between atoms
        totalInteratomPE = 0;
        for (int cellX = 0; cellX < Ncell; ++cellX) {
            for (int cellY = 0; cellY < Ncell; ++cellY) {
                int[] thisCellMembers = cellMembers[cellX][cellY];
                int thisCellPop = cellPops[cellX][cellY];
                for (int off = 0; off < 5; ++off) {
                    int neighborCellX = cellX + cellNeighborOffsetsX[off];
                    int neighborCellY = cellY + cellNeighborOffsetsY[off];
                    if (neighborCellX < 0 || neighborCellY < 0 || neighborCellX >= Ncell || neighborCellY >= Ncell) continue;
                    int[] neighborCellMembers = cellMembers[neighborCellX][neighborCellY];
                    int neighborCellPop = cellPops[neighborCellX][neighborCellY];
                    for (int thisAtom = 0; thisAtom < thisCellPop; ++thisAtom) {
                        for (int neighborAtom = 0; neighborAtom < neighborCellPop; ++neighborAtom) {
                            int i = thisCellMembers[thisAtom];
                            int j = neighborCellMembers[neighborAtom];
                            if (off == 0 && neighborAtom >= thisAtom) break; // avoid double-counting forces within a cell
                            
                            double dx = xs[i] - xs[j];
                            double dy = ys[i] - ys[j];
                            double r2 = dx*dx + dy*dy;
                            if (r2 <= r2cutoff) {
                                // V = 1/r^12 - 2/r^6
                                // F = -dV/dr = 12 (1/r^13 - 1/r^7)
                                // F/r = 12 * 1/r^2 * (1/r^12 - 1/r^6) = 12 (1/r^2) (1/r^6) (1/r^6 - 1)
                                /*double rm2 = 1/r2;
                                double rm6 = rm2*rm2*rm2;
                                double Foverr = 12 * rm6 * (rm6 - 1) * rm2;
                                totalInteratomPE += rm6 * (rm6 - 2);
                                double Fx = Foverr * dx;
                                double Fy = Foverr * dy;
                                axs[i] += Fx;
                                ays[i] += Fy; 
                                axs[j] -= Fx;
                                ays[j] -= Fy;*/
                                
                                // repulsive core
                                double rm2 = 1/r2;
                                double rm6 = rm2*rm2*rm2;
                                double repulsiveFoverR = 12 * rm6 * rm6 * rm2;
                                //System.out.format("%d<->%d  repulsive force, r=%f, F=%f\n", i, j, Math.sqrt(r2), repulsiveFoverR*Math.sqrt(r2));
                                totalInteratomPE += rm6 * rm6;
                                double repulsiveFx = repulsiveFoverR * dx;
                                double repulsiveFy = repulsiveFoverR * dy;
                                axs[i] += repulsiveFx;
                                ays[i] += repulsiveFy; 
                                axs[j] -= repulsiveFx;
                                ays[j] -= repulsiveFy;
                                
                                // dipole-dipole force
                                DoChargePair(i, j, +1, +1);
                                DoChargePair(i, j, +1, -1);
                                DoChargePair(i, j, -1, +1);
                                DoChargePair(i, j, -1, -1);
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
        //for (int i = 0; i < N; ++i) {
        //    System.out.format("acc[%d] = (%f,%f)\n", i, axs[i], ays[i]);
        //}
    }
    
    public void ComputeNewVelocities(double dt) {        
        // compute new velocities
        totalAtomKE = 0;
        totalRotationalKE = 0;
        for (int i = 0; i < N; ++i) {
            vxs[i] = halfvxs[i] + 0.5 * dt * axs[i];
            vys[i] = halfvys[i] + 0.5 * dt * ays[i];
            totalAtomKE += 0.5 * (vxs[i] * vxs[i] + vys[i] * vys[i]);
            omegas[i] = halfomegas[i] + 0.5 * dt * torques[i] / momentOfInertia;
            totalRotationalKE += 0.5 * momentOfInertia * omegas[i] * omegas[i];
        }
        wallV = wallHalfV + 0.5 * dt * wallA;
        if (iter % 1000 == 0 || resampleVelocities) {
            sampleVelocitiesFromMaxwellBoltzmann(targetT);
            resampleVelocities = false;
        }
        wallKE = 0.5 * wallMass * wallV * wallV;
    }
    
    public void EndStep(double dt) {
        iter += 1;
        simTime += dt;        
        avgT = 0.999 * avgT + 0.001 * computeTemperature();
        
        totalE = totalAtomKE + totalRotationalKE + totalInteratomPE + totalWallPE + wallKE + pressureE;
    }
    

    // Maxwell-Boltzmann distribution for each velocity component is gaussian with mean 0
    // and 
    void sampleVelocitiesFromMaxwellBoltzmann(double T)
    {
        double std = Math.sqrt(T);
        double stdOmega = Math.sqrt(T/momentOfInertia);
        for (int i = 0; i < vxs.length; ++i) {
            vxs[i] = std * rand.nextGaussian();
            vys[i] = std * rand.nextGaussian();
            omegas[i] = stdOmega * rand.nextGaussian();
        }
        //double wallStd = Math.sqrt(T / wallMass);
        //wallV = wallStd * rand.nextGaussian();
    }
    
    double computeTemperature() 
    {
        double meanVx2 = 0;
        double meanVy2 = 0;
        for (int i = 0; i < N; ++i) {
            meanVx2 += vxs[i] * vxs[i];
            meanVy2 += vys[i] * vys[i];
        }
        meanVx2 /= N;
        meanVy2 /= N;
        return (meanVx2 + meanVy2)/2;
    }
    
    public void Draw(Graphics g, int width, int height, double scale)
    {
        g.setColor(Color.black);
        g.fillRect(0, 0, width, height);
        
        double diameter = scale;
        if (diameter < 2) diameter = 2;
        double refspeed = 2;
        
        g.setColor(Color.blue);
        for (int i = 0; i < N; ++i) {            
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
            g.fillOval(width/2 + (int)(scale * xs[i]-0.5*diameter), height/2 + (int)(scale * ys[i]-0.5*diameter), (int)diameter, (int)diameter);
            g.setColor(Color.white);
            g.drawLine(width/2 + (int)(scale * xs[i]), height/2 + (int)(scale * ys[i]),
                    width/2 + (int)(scale * xs[i] + 0.5 * diameter * Math.cos(thetas[i])), height/2 + (int)(scale * ys[i] + 0.5 * diameter * Math.sin(thetas[i])));
            g.setColor(Color.black);
            g.drawLine(width/2 + (int)(scale * xs[i]), height/2 + (int)(scale * ys[i]),
                    width/2 + (int)(scale * xs[i] - 0.5 * diameter * Math.cos(thetas[i])), height/2 + (int)(scale * ys[i] - 0.5 * diameter * Math.sin(thetas[i])));
                    
        }
        
        g.setColor(Color.white);
        g.drawRect(width/2 - (int)(scale * wallRadius), height/2 - (int)(scale * wallRadius), (int)(2*scale*wallRadius), (int)(2*scale*wallRadius));
        
        g.setColor(Color.white);
        int textY = -7;
        g.drawString(String.format("targetT = %.4f, temp = %.4f", targetT, avgT), 10, textY += 17);
        g.drawString(String.format("pressure = %.8f", pressure), 10, textY += 17);
        g.drawString(String.format("totalE = %.2f", totalE), 10, textY += 17);
        g.drawString(String.format("totalAtomKE = %.2f", totalAtomKE), 10, textY += 17);
        g.drawString(String.format("totalRotationalKE = %.2f", totalRotationalKE), 10, textY += 17);
        g.drawString(String.format("totalInteratomPE = %.2f", totalInteratomPE), 10, textY += 17);
        g.drawString(String.format("totalWallPE = %.2f", totalWallPE), 10, textY += 17);
        g.drawString(String.format("wallKE = %.2f", wallKE), 10, textY += 17);
        g.drawString(String.format("pressureE = %.2f", pressureE), 10, textY += 17);
        g.drawString(String.format("time = %.2f", simTime), 10, textY += 17);
        g.drawString(String.format("volume = %.2f", 4*wallRadius*wallRadius), 10, textY += 17);
    }
}
