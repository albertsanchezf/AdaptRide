/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package adaptride.entities;

/**
 *
 * @author AlbertSanchez
 */
public class NNDataset {
    
    private float  speed;
    private float  mean_acc_x;
    private float  mean_acc_y;
    private float  mean_acc_z;
    private float  std_acc_x;
    private float  std_acc_y;
    private float  std_acc_z;
    private double sma;
    private float  mean_svm;
    private double entropyX;
    private double entropyY;
    private double entropyZ;
    private int    bike_type;
    private int    phone_location;
    private int    incident_type;
    
    public float getSpeed() {
        return speed;
    }

    public void setSpeed(float speed) {
        this.speed = speed;
    }

    public float getMean_acc_x() {
        return mean_acc_x;
    }

    public void setMean_acc_x(float mean_acc_x) {
        this.mean_acc_x = mean_acc_x;
    }

    public float getMean_acc_y() {
        return mean_acc_y;
    }

    public void setMean_acc_y(float mean_acc_y) {
        this.mean_acc_y = mean_acc_y;
    }

    public float getMean_acc_z() {
        return mean_acc_z;
    }

    public void setMean_acc_z(float mean_acc_z) {
        this.mean_acc_z = mean_acc_z;
    }

    public float getStd_acc_x() {
        return std_acc_x;
    }

    public void setStd_acc_x(float std_acc_x) {
        this.std_acc_x = std_acc_x;
    }

    public float getStd_acc_y() {
        return std_acc_y;
    }

    public void setStd_acc_y(float std_acc_y) {
        this.std_acc_y = std_acc_y;
    }

    public float getStd_acc_z() {
        return std_acc_z;
    }

    public void setStd_acc_z(float std_acc_z) {
        this.std_acc_z = std_acc_z;
    }

    public double getSma() {
        return sma;
    }

    public void setSma(double sma) {
        this.sma = sma;
    }

    public float getMean_svm() {
        return mean_svm;
    }

    public void setMean_svm(float mean_svm) {
        this.mean_svm = mean_svm;
    }

    public double getEntropyX() {
        return entropyX;
    }

    public void setEntropyX(double entropyX) {
        this.entropyX = entropyX;
    }
    
    public double getEntropyY() {
        return entropyY;
    }

    public void setEntropyY(double entropyY) {
        this.entropyY = entropyY;
    }

    public double getEntropyZ() {
        return entropyZ;
    }

    public void setEntropyZ(double entropyZ) {
        this.entropyZ = entropyZ;
    }
    public int getBike_type() {
        return bike_type;
    }

    public void setBike_type(int bike_type) {
        this.bike_type = bike_type;
    }

    public int getPhone_location() {
        return phone_location;
    }

    public void setPhone_location(int phone_location) {
        this.phone_location = phone_location;
    }
    
    public int getIncident_type() {
        return incident_type;
    }

    public void setIncident_type(int incident_type) {
        this.incident_type = incident_type;
    }
    
}
