package adaptride;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jtransforms.fft.DoubleFFT_1D;
import adaptride.entities.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.AbstractList;

/**
 *
 * @author AlbertSanchez
 */
public class AdaptRide {

    private static String ridePath = "/Users/AlbertSanchez/Desktop/Post/Rides/VM2_-1213942287";
    private static String outputPath = "/Users/AlbertSanchez/Desktop/Post/WindowedRides/";
    
    private static int avoidSeconds = 3;
    
    private static int MINNUMBEROFREADINGS = 1;
    
    private static int WINDOWFRAME = 6000; //ms
    private static int WINDOWSPLIT = 3;

    static DoubleFFT_1D fft;
    static double[] xFFT;
    static double[] yFFT;
    static double[] zFFT;
    static double[] fftDataX;
    static double[] fftDataY;
    static double[] fftDataZ;

    public static void main(String[] args) 
    {
        
        // Ensure final '/' in path strings
        if(!outputPath.endsWith("/")) outputPath += "/";
        
        try
        {
            Ride r = avoidInitAndEndSeconds(getRide(ridePath),avoidSeconds);
            
            if (r != null)
            {
                List<WindowedRide> wr = chopRide(r);
                List<NNDataset> ds = calculateStatistics(wr);
                writeCSVFile(outputPath, ds);
            }
        }
        catch (Exception e)
        {
            System.out.println(e.getMessage());
        }

    }
    
    public static Ride getRide(String ridePath) throws IOException
    {
        boolean fileWithWrongFormat = false;
        Ride ride;

        String[] incidentFields;
                
        List<Double> latitude, longitude, 
                     tmp_latitude, tmp_longitude;
        List<Float> acc_x, acc_y,  acc_z, acc_68, gyr_a, gyr_b, gyr_c,
                    tmp_acc_x, tmp_acc_y, tmp_acc_z, tmp_acc_68, tmp_gyr_a, tmp_gyr_b, tmp_gyr_c;
        List<Long> timestamp, tmp_timestamp;

        double prevLat = 0, prevLon = 0;
        float  prevAcc_68 = 0, prevGyr_a = 0, prevGyr_b = 0, prevGyr_c = 0;
        int    readings = 0;

        int l_number = 0;
        List<Integer> tmp_linenumber = new ArrayList<>();
        List<Integer> linenumber = new ArrayList<>();

        fileWithWrongFormat = false;
        ride = new Ride();
        ride.setDs_name(ridePath);

        latitude = new ArrayList<>();
        longitude = new ArrayList<>();
        acc_x = new ArrayList<>();
        acc_y = new ArrayList<>();
        acc_z = new ArrayList<>();
        timestamp = new ArrayList<>();
        acc_68 = new ArrayList<>();
        gyr_a = new ArrayList<>();
        gyr_b = new ArrayList<>();
        gyr_c = new ArrayList<>();

        tmp_latitude = new ArrayList<>();
        tmp_longitude = new ArrayList<>();
        tmp_acc_x = new ArrayList<>();
        tmp_acc_y = new ArrayList<>();
        tmp_acc_z = new ArrayList<>();
        tmp_timestamp = new ArrayList<>();
        tmp_acc_68 = new ArrayList<>();
        tmp_gyr_a = new ArrayList<>();
        tmp_gyr_b = new ArrayList<>();
        tmp_gyr_c = new ArrayList<>();

        FileReader r = new FileReader(ridePath);
        BufferedReader br = new BufferedReader(r);
        
        br.readLine(); //Reading version (e.g. <i5#1> in iOS or 45#2 in android)
        br.readLine(); //Reading headers (key,lat,lon,ts,bike,childCheckBox,...)
        
        l_number+=2;
        
        String line = br.readLine();
        l_number++;
        if(line.startsWith("0"))
        {
            incidentFields = line.split(",",-1); 

            ride.setBikeType(Integer.parseInt(incidentFields[4]));
            ride.setPhoneLocation(Integer.parseInt(incidentFields[7]));

            // Read until header of the ride data
            line = br.readLine();
            l_number++;
            while(!line.equals("lat,lon,X,Y,Z,timeStamp,acc,a,b,c"))
            {
                line = br.readLine();
                l_number++;
                if(line.equals("lat,lon,X,Y,Z,timeStamp")) 
                {
                    fileWithWrongFormat = true;
                    break;
                }
            }

            if (!fileWithWrongFormat) 
            {
                line = br.readLine();
                l_number++;
                if (line != null)
                {
                    incidentFields = line.split(",",-1); 

                    // Save values that are not in constant updating
                    prevLat = Double.parseDouble(incidentFields[0]);
                    prevLon = Double.parseDouble(incidentFields[1]);
                    prevAcc_68 = Float.parseFloat(incidentFields[6]);
                    prevGyr_a  = Float.parseFloat(incidentFields[7]);
                    prevGyr_b  = Float.parseFloat(incidentFields[8]);
                    prevGyr_c  = Float.parseFloat(incidentFields[9]); 

                    // Add the first line values
                    tmp_linenumber.add(l_number);
                    tmp_latitude.add(prevLat);
                    tmp_longitude.add(prevLon);
                    tmp_acc_x.add(Float.parseFloat(incidentFields[2]));
                    tmp_acc_y.add(Float.parseFloat(incidentFields[3]));
                    tmp_acc_z.add(Float.parseFloat(incidentFields[4]));
                    tmp_timestamp.add(Long.parseLong(incidentFields[5]));
                    tmp_acc_68.add(prevAcc_68);
                    tmp_gyr_a.add(prevGyr_a);
                    tmp_gyr_b.add(prevGyr_b);
                    tmp_gyr_c.add(prevGyr_c);
                    line = br.readLine();
                    l_number++;

                    // Read all the file
                    while (line != null)
                    {
                        incidentFields = line.split(",",-1); 
                        // Counter for readings
                        // If there is no data in Latitude field (incidentFields[0]) we
                        // add the counter
                        if (incidentFields[0].equals(""))
                        {
                            readings++;
                        }
                        // If there is data we save the tmp data if the # of readings is enough
                        else
                        {
                            // If the number of readings is sufficient we save it
                            if (readings >= MINNUMBEROFREADINGS)
                            {
                                linenumber.addAll(tmp_linenumber);
                                latitude.addAll(tmp_latitude);
                                longitude.addAll(tmp_longitude);
                                acc_x.addAll(tmp_acc_x);
                                acc_y.addAll(tmp_acc_y);
                                acc_z.addAll(tmp_acc_z);
                                timestamp.addAll(tmp_timestamp);
                                acc_68.addAll(tmp_acc_68);
                                gyr_a.addAll(tmp_gyr_a);
                                gyr_b.addAll(tmp_gyr_b);
                                gyr_c.addAll(tmp_gyr_c);
                            }
                            tmp_linenumber.clear();
                            tmp_latitude.clear();
                            tmp_longitude.clear();
                            tmp_acc_x.clear();
                            tmp_acc_y.clear();
                            tmp_acc_z.clear();
                            tmp_timestamp.clear();
                            tmp_acc_68.clear();
                            tmp_gyr_a.clear();
                            tmp_gyr_b.clear();
                            tmp_gyr_c.clear();

                            readings = 1;
                        }

                        tmp_linenumber.add(l_number);
                        
                        if (incidentFields[0].equals(""))
                            tmp_latitude.add(prevLat);
                        else
                        {
                            tmp_latitude.add(Double.parseDouble(incidentFields[0]));
                            prevLat = tmp_latitude.get(tmp_latitude.size()-1);
                        }

                        if (incidentFields[1].equals(""))
                            tmp_longitude.add(prevLon);
                        else
                        {
                            tmp_longitude.add(Double.parseDouble(incidentFields[1]));
                            prevLon = tmp_longitude.get(tmp_longitude.size()-1);
                        }

                        tmp_acc_x.add(Float.parseFloat(incidentFields[2]));
                        tmp_acc_y.add(Float.parseFloat(incidentFields[3]));
                        tmp_acc_z.add(Float.parseFloat(incidentFields[4]));
                        tmp_timestamp.add(Long.parseLong(incidentFields[5]));
                        if (incidentFields[6].equals(""))
                            tmp_acc_68.add(prevAcc_68);
                        else
                        {
                            tmp_acc_68.add(Float.parseFloat(incidentFields[6]));
                            prevAcc_68 = tmp_acc_68.get(tmp_acc_68.size()-1);
                        }
                        if (incidentFields[7].equals(""))
                            tmp_gyr_a.add(prevGyr_a);
                        else
                        {
                            tmp_gyr_a.add(Float.parseFloat(incidentFields[7]));
                            prevGyr_a = tmp_gyr_a.get(tmp_gyr_a.size()-1);
                        }
                        if (incidentFields[8].equals(""))
                            tmp_gyr_b.add(prevGyr_b);
                        else
                        {
                            tmp_gyr_b.add(Float.parseFloat(incidentFields[8]));
                            prevGyr_b = tmp_gyr_b.get(tmp_gyr_b.size()-1);
                        }
                        if (incidentFields[9].equals(""))
                            tmp_gyr_c.add(prevGyr_c);
                        else
                        {
                            tmp_gyr_c.add(Float.parseFloat(incidentFields[9]));
                            prevGyr_c = tmp_gyr_c.get(tmp_gyr_c.size()-1);
                        }

                        // Read next line
                        line = br.readLine();
                        l_number++;

                        if (line == null && readings >= MINNUMBEROFREADINGS)
                        {
                            linenumber.addAll(tmp_linenumber);
                            tmp_linenumber.clear();
                            latitude.addAll(tmp_latitude);
                            longitude.addAll(tmp_longitude);
                            acc_x.addAll(tmp_acc_x);
                            acc_y.addAll(tmp_acc_y);
                            acc_z.addAll(tmp_acc_z);
                            timestamp.addAll(tmp_timestamp);
                            acc_68.addAll(tmp_acc_68);
                            gyr_a.addAll(tmp_gyr_a);
                            gyr_b.addAll(tmp_gyr_b);
                            gyr_c.addAll(tmp_gyr_c);
                            tmp_latitude.clear();
                            tmp_longitude.clear();
                            tmp_acc_x.clear();
                            tmp_acc_y.clear();
                            tmp_acc_z.clear();
                            tmp_timestamp.clear();
                            tmp_acc_68.clear();
                            tmp_gyr_a.clear();
                            tmp_gyr_b.clear();
                            tmp_gyr_c.clear();
                            readings = 1;
                        }
                    }

                    // Saving detailed data to the ride
                    ride.setLinenumber(linenumber);
                    
                    if(!latitude.isEmpty())
                        ride.setLatitude(latitude);
                    if(!longitude.isEmpty())
                        ride.setLongitude(longitude);
                    if(!acc_x.isEmpty())
                        ride.setAcc_x(acc_x);
                    if(!acc_y.isEmpty())
                        ride.setAcc_y(acc_y);
                    if(!acc_z.isEmpty())
                        ride.setAcc_z(acc_z);
                    if(!timestamp.isEmpty())
                        ride.setTimestamp(timestamp);
                    if(!acc_68.isEmpty())
                        ride.setAcc_68(acc_68);
                    if(!gyr_a.isEmpty())
                        ride.setGyr_a(gyr_a);
                    if(!gyr_b.isEmpty())
                        ride.setGyr_b(gyr_b);
                    if(!gyr_c.isEmpty())
                        ride.setGyr_c(gyr_c);

                }
                else 
                {
                    System.out.println("Filename: " + ridePath + " is empty");
                    ride = null;
                }
            }
            else 
            {
                System.out.println("Filename: " + ridePath + " has a wrong format");
                System.out.println("Missing ride data (acc_68,gyr_a,gyr_b,gyr_c");
                ride = null;
            }
        }
        else
        {
            System.out.println("Filename: " + ridePath + " has a wrong format");
            System.out.println("Neither bike type nor phone location set");
            ride = null;
        }
        return ride;

    }
    
    public static List<WindowedRide> chopRide(Ride r)
    {
        List<WindowedRide> windowedRides = new ArrayList<>();
        WindowedRide ride;
        int i=0;
        long[] timestamps;
        List<Integer> initIndexes, endIndexes;
        Long dt = 0l;
        int bikeType = 0, phoneLocation = 0;
          
        System.out.println(r.getDs_name());
        initIndexes = new ArrayList<>();
        endIndexes = new ArrayList<>();

        // Phone Location and Bike Type should not change between incidents
        bikeType = r.getBikeType();
        phoneLocation = r.getPhoneLocation();
        
        // We check that the ride is longer than the actual Windowframe span
        if(r.getTimestamp().get(r.getTimestamp().size()-1) - r.getTimestamp().get(0) >= WINDOWFRAME)
        {
            timestamps = r.getTimestamp().stream().mapToLong(x->x).toArray();
            initIndexes.add(i);

            for (int j=1; j < timestamps.length; j++)
            {
                dt = timestamps[j] - timestamps[i];
                if (dt >= WINDOWFRAME)
                {
                    endIndexes.add(j);
                    j=i;
                    while((timestamps[j]-timestamps[i])< WINDOWFRAME/WINDOWSPLIT)
                    {
                        j++;
                        System.out.println(j);
                        if (j==timestamps.length-1)
                            break; 
                   }
                    j--;
                    i=j;
                    initIndexes.add(i);
                }
                else if (j==timestamps.length-1)
                    endIndexes.add(j);
                
            }
            int ii = 0, ei = 0;
            for (int k=0; k < initIndexes.size(); k++)
            {
                ii = initIndexes.get(k);
                ei = endIndexes.get(k);

                ride = new WindowedRide();

                ride.setIncident(0); //Not relevant because we want to predict                

                ride.setDs_name(r.getDs_name());
                ride.setPhoneLocation(phoneLocation);
                ride.setBikeType(bikeType);
                ride.setLatitude(r.getLatitude().subList(ii, ei));
                ride.setLongitude(r.getLongitude().subList(ii, ei));
                ride.setAcc_x(r.getAcc_x().subList(ii, ei));
                ride.setAcc_y(r.getAcc_y().subList(ii, ei));
                ride.setAcc_z(r.getAcc_z().subList(ii, ei));
                ride.setTimestamp(r.getTimestamp().subList(ii, ei));
                ride.setAcc_68(r.getAcc_68().subList(ii, ei));
                ride.setGyr_a(r.getGyr_a().subList(ii, ei));
                ride.setGyr_b(r.getGyr_b().subList(ii, ei));
                ride.setGyr_c(r.getGyr_c().subList(ii, ei));

                windowedRides.add(ride);
            }
        }
        return windowedRides;

    }
    
    public static List<NNDataset> calculateStatistics(List<WindowedRide> rides)
    {
        List<NNDataset> nnDatasetList = new ArrayList<>();
        NNDataset nnDataset;
        long t0 = 0l, t1 = 0l, d = 0l;
        int numberOfNNDataset = 0;
        
        for (WindowedRide ride : rides)
        {
            nnDataset = new NNDataset();
            
            nnDataset.setSpeed(getSpeed(ride));
            nnDataset.setMean_acc_x(getMeanF(ride.getAcc_x()));
            nnDataset.setMean_acc_y(getMeanF(ride.getAcc_y()));
            nnDataset.setMean_acc_z(getMeanF(ride.getAcc_z()));
            nnDataset.setStd_acc_x(getStdDeviation(ride.getAcc_x()));
            nnDataset.setStd_acc_y(getStdDeviation(ride.getAcc_y()));
            nnDataset.setStd_acc_z(getStdDeviation(ride.getAcc_z()));
            nnDataset.setSma(getSma(ride.getAcc_x(),ride.getAcc_y(),ride.getAcc_z()));
            nnDataset.setMean_svm(getMeanD(getSvm(ride.getAcc_x(),ride.getAcc_y(),ride.getAcc_z())));
                        
            calculateFFT(ride);
            nnDataset.setEntropyX(getEntropy(xFFT));
            nnDataset.setEntropyY(getEntropy(yFFT));
            nnDataset.setEntropyZ(getEntropy(zFFT));
                        
            nnDataset.setBike_type(ride.getBikeType());
            nnDataset.setPhone_location(ride.getPhoneLocation());
            nnDataset.setIncident_type(ride.getIncident());
            
            nnDatasetList.add(nnDataset);
            numberOfNNDataset++;
        }
        
        System.out.println(numberOfNNDataset + " NNDataset records added");
        return nnDatasetList;
    }
    
    public static float getSpeed(WindowedRide ride) // m/s
    {
        double iniLat = ride.getLatitude().get(0);
        double iniLon = ride.getLongitude().get(0);
        double endLat = ride.getLatitude().get(ride.getLatitude().size()-1);
        double endLon = ride.getLongitude().get(ride.getLongitude().size()-1);
        long dt = (ride.getTimestamp().get(ride.getTimestamp().size()-1) - ride.getTimestamp().get(0))/1000; // s
        
        double earthRadius = 6371000; //m (3958.75 miles)
        double dLat = Math.toRadians(endLat - iniLat);
        double dLon = Math.toRadians(endLon - iniLon);
        double sindLat = Math.sin(dLat/2);
        double sindLon = Math.sin(dLon/2);
        double a = Math.pow(sindLat,2) + Math.pow(sindLon,2)
                * Math.cos(Math.toRadians(iniLat)) * Math.cos(Math.toRadians(endLat));
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
        
        return (float) ((earthRadius * c) / dt);
    }
    
    public static float getMeanF(List<Float> x)
    {
        return (float) x.stream().mapToDouble(i->i).average().orElse(0);
    }
    
    public static float getMeanD(List<Double> x)
    {
        return (float) x.stream().mapToDouble(i->i).average().orElse(0);
    }
    
    public static double getSma(List<Float> x, List<Float> y, List<Float> z)
    {
        double sma = 0;
        
        double xA[] = x.stream().mapToDouble(i->i).toArray();
        double yA[] = y.stream().mapToDouble(i->i).toArray();
        double zA[] = z.stream().mapToDouble(i->i).toArray();

        if (xA.length == yA.length && yA.length == zA.length)
        {
            for (int i = 0; i < xA.length; i++)
            {
                sma += Math.abs(xA[i]) + Math.abs(yA[i]) + Math.abs(zA[i]);
            }
        }

        return sma/xA.length;
    }
    
    public static List<Double> getSvm(List<Float> x, List<Float> y, List<Float> z)
    {
        int i=0;
        List<Double> svm = new ArrayList<>();
        
        for (Float f : x)
        {
            double d = Math.sqrt(Math.pow(f,2)) + Math.sqrt(Math.pow(y.get(i),2)) + Math.sqrt(Math.pow(z.get(i),2));
            svm.add(d);
            i++;
        }
        
        return svm;
        
    }
    
    public static float getMean_svm(List<Float> x, List<Float> y, List<Float> z)
    {
        
        double xA[] = x.stream().mapToDouble(i->i).toArray();
        double yA[] = y.stream().mapToDouble(i->i).toArray();
        double zA[] = z.stream().mapToDouble(i->i).toArray();
        
        double svm[] = new double[x.size()];
        
        for (int i=0; i < xA.length; i++)
        {
            svm[i] = Math.sqrt(Math.pow(xA[i], 2) + Math.pow(yA[i], 2) + Math.pow(zA[i], 2));
        }

        return (float) Arrays.stream(svm).average().orElse(0);
    }
    
    public static double getVariance(List<Float> x)
    {
        double variance = 0;
        
        double xA[] = x.stream().mapToDouble(i->i).toArray();

        if (xA.length > 0)
        {
            float mean = getMeanF(x);

            for (int i = 0; i < xA.length; i++)
            {
                    variance += Math.pow(xA[i] - mean, 2);
            }

            variance = variance / (float) xA.length;
        }

        return variance;
    }

    public static float getStdDeviation(List<Float> x)
    {
        float stdDeviation = 0;

        if (x.size() > 0)
        {
            stdDeviation = (float) Math.sqrt(getVariance(x));
        }

        return stdDeviation;
    }
    
    public static double getEnergy(double[] fftValues)
    {
        float energy = 0;

        // Calculate energy
        for (int i = 0; i < fftValues.length; i++)
        {
            energy += Math.pow(fftValues[i], 2);
        }

        if (fftValues.length > 0)
        {
            energy = energy / (float) fftValues.length;
        }

        return energy;
    }

    private static double getEntropy(double[] fftValues)
    {
        double entropy = 0;
        double[] psd = new double[fftValues.length];

        if (fftValues.length > 0)
        {
            // Calculate Power Spectral Density
            for (int i = 0; i < fftValues.length; i++)
            {
                psd[i] = (Math.pow(fftValues[i], 2) / fftValues.length);
            }

            double psdSum = getSum(psd);

            if (psdSum > 0)
            {
                // Normalize calculated PSD so that it can be viewed as a Probability Density Function
                for (int i = 0; i < fftValues.length; i++)
                {
                    psd[i] = psd[i] / psdSum;
                }

                // Calculate the Frequency Domain Entropy
                for (int i = 0; i < fftValues.length; i++)
                {
                    if (psd[i] != 0)
                    {
                        entropy += psd[i] * Math.log(psd[i]);
                    }
                }

                entropy *= -1;
            }
        }

    return entropy;
    
    }
    
    public static double getSum(double[] values)
    {
        double sum = 0;

        for (int i = 0; i < values.length; i++)
        {
            sum += values[i];
        }

        return sum;
    }
    
    public static void calculateFFT(WindowedRide r)
    {
        xFFT = r.getAcc_x().stream().mapToDouble(i->i).toArray();
        yFFT = r.getAcc_y().stream().mapToDouble(i->i).toArray();
        zFFT = r.getAcc_z().stream().mapToDouble(i->i).toArray();
        
        fft = new DoubleFFT_1D(xFFT.length);
        
        fftDataX = new double[2*xFFT.length];
        fftDataY = new double[2*yFFT.length];
        fftDataZ = new double[2*zFFT.length];
        
        for (int i = 0; i < xFFT.length; i++) 
        {
            // copying data to the fft data buffer, imaginary part is 0
            fftDataX[2 * i] = xFFT[i];
            fftDataX[2 * i + 1] = 0;
            fftDataY[2 * i] = yFFT[i];
            fftDataY[2 * i + 1] = 0;
            fftDataZ[2 * i] = zFFT[i];
            fftDataZ[2 * i + 1] = 0;
        }
        
        fft.complexForward(fftDataX);
        fft.complexForward(fftDataY);
        fft.complexForward(fftDataZ);
        
    }
    
    public static String writeCSVFile(String path, List<NNDataset> nnDataset) throws IOException
    {
        String filename = path + "ride.csv";
        String line = "";
        String[] s, s1;
        int csvRecords = 0, maxFile=0;
        boolean dsFile = false;
        
        // To ensure any file is overwrited
        File f = new File(path);
        String[] files = f.list();
        for(int i=0; i < files.length; i++)
        {
            if(new File(path + files[i]).isFile())
            {
                s = files[i].split("ride"); 
                if(s.length==2)
                {
                    if(s[1].equals(".csv"))
                        dsFile = true;
                    else
                    {
                        s1 = s[1].split(".csv");
                        if (s1[0].matches("^[0-9]*$"))
                            if(Integer.valueOf(s1[0]) > maxFile)
                                maxFile = Integer.valueOf(s1[0]);            
                    }
                } 
            }
        }
        
        maxFile++;
        
        if (dsFile)
            new File(filename).renameTo(new File(path + "ride" + maxFile + ".csv"));
                
        FileWriter writer = new FileWriter(filename);
        
        //line = "speed,mean_acc_x,mean_acc_y,mean_acc_z,std_acc_x,std_acc_y,std_acc_z,sma,mean_svm,entropyX,entropyY,entropyZ,bike_type,phone_location,incident_type\n";
        //writer.append(line);
        
        for(NNDataset l : nnDataset)
        {
            line  = l.getSpeed() + ",";
            line += l.getMean_acc_x() + ",";
            line += l.getMean_acc_y() + ",";
            line += l.getMean_acc_z() + ",";
            line += l.getStd_acc_x() + ",";
            line += l.getStd_acc_y() + ",";
            line += l.getStd_acc_z() + ",";
            line += l.getSma() + ",";
            line += l.getMean_svm() + ",";
            line += l.getEntropyX() + ",";
            line += l.getEntropyY() + ",";
            line += l.getEntropyZ() + ",";
            line += l.getBike_type()  + ",";
            line += l.getPhone_location() + ",";
            line += l.getIncident_type() + "\n";
            writer.append(line);            
            csvRecords++;   
        }
        
        writer.flush();
        writer.close();
        System.out.println(String.format("Csv records: %d",csvRecords));
        return filename;
    }
    
    public static Ride avoidInitAndEndSeconds(Ride ride, int seconds)
    {
        Ride r = new Ride();
        long t0 = 0l, t1 = 0l;
        long[] timestamps;
        int initIndex = 0, endIndex = 0;

        r.setDs_name(ride.getDs_name());
        r.setPhoneLocation(ride.getPhoneLocation());
        r.setBikeType(ride.getBikeType());
        
        timestamps = ride.getTimestamp().stream().mapToLong(x->x).toArray();
        endIndex = timestamps.length-1;
        
        t0 = timestamps[0] + (seconds * 1000);
        t1 = timestamps[timestamps.length-1] - (seconds * 1000);
        
        for (int i=0; i<ride.getLinenumber().size(); i++)
        {
            if (timestamps[i] < t0)
                initIndex = i;
            
            if(timestamps[i] < t1)
                endIndex = i;
        }
        
        initIndex++;
        
        r.setLinenumber(ride.getLinenumber().subList(initIndex, endIndex));
        r.setLatitude(ride.getLatitude().subList(initIndex, endIndex));
        r.setLongitude(ride.getLongitude().subList(initIndex, endIndex));
        r.setAcc_x(ride.getAcc_x().subList(initIndex, endIndex));
        r.setAcc_y(ride.getAcc_y().subList(initIndex, endIndex));
        r.setAcc_z(ride.getAcc_z().subList(initIndex, endIndex));
        r.setAcc_68(ride.getAcc_68().subList(initIndex, endIndex));
        r.setGyr_a(ride.getGyr_a().subList(initIndex, endIndex));
        r.setGyr_b(ride.getGyr_b().subList(initIndex, endIndex));
        r.setGyr_c(ride.getGyr_c().subList(initIndex, endIndex));
        r.setTimestamp(ride.getTimestamp().subList(initIndex, endIndex));
        
        return r;
    }
}
