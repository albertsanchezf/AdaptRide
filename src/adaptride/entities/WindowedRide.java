/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package adaptride.entities;

import java.util.List;

/**
 *
 * @author AlbertSanchez
 */
public class WindowedRide {
    
    private String ds_name;
    private int bikeType;
    private int phoneLocation;
    private int incident;
    private List<Double> latitude;
    private List<Double> longitude;
    private List<Float> acc_x;
    private List<Float> acc_y;
    private List<Float> acc_z;
    private List<Long> timestamp;
    private List<Float> acc_68;
    private List<Float> gyr_a;
    private List<Float> gyr_b;
    private List<Float> gyr_c;

    public String getDs_name() {
        return ds_name;
    }

    public void setDs_name(String ds_name) {
        this.ds_name = ds_name;
    }

    public int getBikeType() {
        return bikeType;
    }

    public void setBikeType(int bikeType) {
        this.bikeType = bikeType;
    }    
    
    public int getPhoneLocation() {
        return phoneLocation;
    }

    public void setPhoneLocation(int phoneLocation) {
        this.phoneLocation = phoneLocation;
    }        
    
    public int getIncident() {
        return incident;
    }

    public void setIncident(int incident) {
        this.incident = incident;
    }    
    
    public List<Double> getLatitude() {
        return latitude;
    }

    public void setLatitude(List<Double> latitude) {
        this.latitude = latitude;
    }

    public List<Double> getLongitude() {
        return longitude;
    }

    public void setLongitude(List<Double> longitude) {
        this.longitude = longitude;
    }

    public List<Float> getAcc_x() {
        return acc_x;
    }

    public void setAcc_x(List<Float> acc_x) {
        this.acc_x = acc_x;
    }

    public List<Float> getAcc_y() {
        return acc_y;
    }

    public void setAcc_y(List<Float> acc_y) {
        this.acc_y = acc_y;
    }

    public List<Float> getAcc_z() {
        return acc_z;
    }

    public void setAcc_z(List<Float> acc_z) {
        this.acc_z = acc_z;
    }

    public List<Long> getTimestamp() {
        return timestamp;
    }

    public void setTimestamp(List<Long> timestamp) {
        this.timestamp = timestamp;
    }

    public List<Float> getAcc_68() {
        return acc_68;
    }

    public void setAcc_68(List<Float> acc_68) {
        this.acc_68 = acc_68;
    }

    public List<Float> getGyr_a() {
        return gyr_a;
    }

    public void setGyr_a(List<Float> gyr_a) {
        this.gyr_a = gyr_a;
    }

    public List<Float> getGyr_b() {
        return gyr_b;
    }

    public void setGyr_b(List<Float> gyr_b) {
        this.gyr_b = gyr_b;
    }

    public List<Float> getGyr_c() {
        return gyr_c;
    }

    public void setGyr_c(List<Float> gyr_c) {
        this.gyr_c = gyr_c;
    }
    
}
