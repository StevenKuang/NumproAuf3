import java.util.Arrays;

public class Gauss {

    /**
     * Diese Methode soll die Loesung x des LGS R*x=b durch
     * Rueckwaertssubstitution ermitteln.
     * PARAMETER:
     * R: Eine obere Dreiecksmatrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] backSubst(double[][] R, double[] b) {
        //TODO: Diese Methode ist zu implementieren
        return null;
    }

    /**
     * Diese Methode soll die Loesung x des LGS A*x=b durch Gauss-Elimination mit
     * Spaltenpivotisierung ermitteln. A und b sollen dabei nicht veraendert werden.
     * PARAMETER: A:
     * Eine regulaere Matrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     *///
    public static double[] solve(double[][] A, double[] b) {
        double[][] temp = Arrays.stream(A).map(double[]::clone).toArray(double[][]::new);

        double[] temp1 = Arrays.copyOf(b, b.length);;
        int length = b.length;
        for (int pivot = 0; pivot < length; pivot++) {
            int max = pivot;
            for (int i = pivot + 1; i < length; i++) {
                if (Math.abs(temp[i][pivot]) > Math.abs(temp[max][pivot])) {
                    max = i;
                }
            }
            double swap2 = temp1[pivot]; temp1[pivot] = temp1[max]; temp1[max] = swap2;
            double[] swap = temp[pivot]; temp[pivot] = temp[max]; temp[max] = swap;
            // Gauss
            for (int i = pivot + 1; i < length; i++) {
                double Coefficient = temp[i][pivot] / temp[pivot][pivot];
                temp1[i] -= Coefficient * temp1[pivot];
                for (int j = pivot; j < length; j++) {
                    temp[i][j] -= Coefficient * temp[pivot][j];
                }
            }
        }
        //solve equations
        double[] x = new double[length];
        for (int i = length - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < length; j++) {
                sum += temp[i][j] * x[j];
            }
            x[i] = (temp1[i] - sum) / temp[i][i];
        }

        return x;
    }

    /**
     * Diese Methode soll eine Loesung p!=0 des LGS A*p=0 ermitteln. A ist dabei
     * eine nicht invertierbare Matrix. A soll dabei nicht veraendert werden.
     *
     * Gehen Sie dazu folgendermassen vor (vgl.Aufgabenblatt):
     * -Fuehren Sie zunaechst den Gauss-Algorithmus mit Spaltenpivotisierung
     *  solange durch, bis in einem Schritt alle moeglichen Pivotelemente
     *  numerisch gleich 0 sind (d.h. <1E-10)
     * -Betrachten Sie die bis jetzt entstandene obere Dreiecksmatrix T und
     *  loesen Sie Tx = -v durch Rueckwaertssubstitution
     * -Geben Sie den Vektor (x,1,0,...,0) zurueck
     *
     * Sollte A doch intvertierbar sein, kann immer ein Pivot-Element gefunden werden(>=1E-10).
     * In diesem Fall soll der 0-Vektor zurueckgegeben werden.
     * PARAMETER:
     * A: Eine singulaere Matrix der Groesse n x n
     */
    public static double[] solveSing(double[][] A) {
        //TODO: Diese Methode ist zu implementieren
        return null;
    }

    /**
     * Diese Methode berechnet das Matrix-Vektor-Produkt A*x mit A einer nxm
     * Matrix und x einem Vektor der Laenge m. Sie eignet sich zum Testen der
     * Gauss-Loesung
     */
    public static double[] matrixVectorMult(double[][] A, double[] x) {
        int n = A.length;
        int m = x.length;

        double[] y = new double[n];

        for (int i = 0; i < n; i++) {
            y[i] = 0;
            for (int j = 0; j < m; j++) {
                y[i] += A[i][j] * x[j];
            }
        }

        return y;
    }
}
