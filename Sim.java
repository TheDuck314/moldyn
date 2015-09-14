import java.awt.Graphics;


public interface Sim {
    public void ComputeHalfVelocities(double dt);
    public void ComputeNewPositions(double dt);
    public void ComputeAcclerations(double dt);
    public void ComputeNewVelocities(double dt);
    public void EndStep(double dt);
    public void Draw(Graphics g, int width, int height, double scale);
    
    public void SetT(double T);
    public void SetP(double P);
    public void ResampleVelocities();
}
