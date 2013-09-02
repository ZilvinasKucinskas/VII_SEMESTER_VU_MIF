import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * @author Žilvinas Kučinskas
 *
 * Aprasymas:
 *
 * Klase skirta rasyti i faila binarini vektoriu kaip stringa
 *
 */

public class BinaryMatrixWriter {
    private File file;

    /*
     * Konstruktorius, kuris kaip parametra gauna faila i kuri rasysime
     */
    public BinaryMatrixWriter(File file) {
        this.file = file;
    }

    /*
     * Metodas skirtas rasyti i duotaji faila vektoriu
     */
    public boolean writeVector(int[][] vector) {
        try {
        PrintWriter pw = new PrintWriter(new FileWriter(file));
        for (int i = 0; i < vector[0].length; i++) {
            pw.print(vector[0][i]);
        }
        
        pw.close();
        } catch (Exception e) {
            System.out.println("Error: " + e);
            return false;
        }
        return true;
    }
}