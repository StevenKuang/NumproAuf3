import java.util.Arrays;

public class Gauss {

    public static final double EPSILON = 1E-10;
    /**
     * Diese Methode soll die Loesung x des LGS R*x=b durch
     * Rueckwaertssubstitution ermitteln.
     * PARAMETER:
     * R: Eine obere Dreiecksmatrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] backSubst(double[][] R, double[] b) {
        double[] result = new double[R.length];
        for (int i = R.length - 1; i >= 0; i--) {
            double tempValue = 0.0;
            for (int j = i + 1; j < R.length; j++) {
                tempValue += R[i][j] * result[j];
            }
            result[i] = (b[i] - tempValue) / R[i][i];
        }

        return result;
    }

    /**
     * Diese Methode soll die Loesung x des LGS A*x=b durch Gauss-Elimination mit
     * Spaltenpivotisierung ermitteln. A und b sollen dabei nicht veraendert werden.
     * PARAMETER: A:
     * Eine regulaere Matrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     *///
    public static double[] solve(double[][] A, double[] b) {
        double[][] tempA = Arrays.stream(A).map(double[]::clone).toArray(double[][]::new);

        double[] tempB = Arrays.copyOf(b, b.length);
        int length = b.length;
        for (int pivot = 0; pivot < length; pivot++) {
            int maxPivot = pivot;
            for (int i = pivot + 1; i < length; i++) {
                if (Math.abs(tempA[i][pivot]) > Math.abs(tempA[maxPivot][pivot])) {
                    maxPivot = i;
                }
            }
            double swapVector = tempB[pivot]; tempB[pivot] = tempB[maxPivot]; tempB[maxPivot] = swapVector;
            double[] swapMatrix = tempA[pivot]; tempA[pivot] = tempA[maxPivot]; tempA[maxPivot] = swapMatrix;
            // Gauss
            for (int i = pivot + 1; i < length; i++) {
                double Coefficient = tempA[i][pivot] / tempA[pivot][pivot];
                tempB[i] -= Coefficient * tempB[pivot];
                for (int j = pivot; j < length; j++) {
                    tempA[i][j] -= Coefficient * tempA[pivot][j];
                }
            }
        }
        return backSubst(tempA, tempB);
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
        double[][] copyA = Arrays.stream(A).map(double[]::clone).toArray(double[][]::new);
        boolean noPivot = false;
        int length = A.length; int pivot;
        for (pivot = 0; pivot < length; pivot++) {
            int maxPivot = pivot;
            for (int i = pivot + 1; i < length; i++) {
                if (Math.abs(copyA[i][pivot]) > Math.abs(copyA[maxPivot][pivot])) {
                    maxPivot = i;
                }
            }
            if (copyA[maxPivot][pivot] < EPSILON) { noPivot = true; break;}
            double[] swapMatrix = copyA[pivot]; copyA[pivot] = copyA[maxPivot]; copyA[maxPivot] = swapMatrix;
            // Gauss
            for (int i = pivot + 1; i < length; i++) {
                double Coefficient = copyA[i][pivot] / copyA[pivot][pivot];
                for (int j = pivot; j < length; j++) {
                    copyA[i][j] -= Coefficient * copyA[pivot][j];
                }
            }
        }
        if(noPivot) {
            if(copyA.length == 2) return new double[] {-copyA[0][1]/copyA[0][0], 1.0};
            // extract the vector v which the column where we stopped the gauss elimination => pivot-1
            double[] v = new double[pivot];
            for (int i = 0; i <pivot  ; i++) {
                v[i] = - copyA[i][pivot];
            }
            // extract the matrix which is the part of the original matrix from 0->pivot-1
            double[][] T = new double[pivot][pivot];
            for (int i = 0; i < pivot; i++) {
                System.arraycopy(copyA[i], 0, T[i], 0, pivot);
            }
            double[] x1 = backSubst(T,v);

            double[] x = new double[A.length];

            // the last index that was filled was k
            for (int i = 0; i < x1.length; i++) {
                x[i] = x1[i];
            }
            x[pivot] = 1.0;
            for (int i = pivot+1; i < x.length; i++) {
                x[i] = 0;
            }
            return x;
        }
        //if a is singular, which means we always found a pivot , we return the null vector

        return new double[A.length];
    }
    private static int getIndexRowPivotUnderElement(double[][] A, int column, int row) {
        //A[0].length is number of rows, and A.length is number of columns
        int max = row; //max is the index of the row that contains the pivot
        for (int i = row+1; i < A.length ; i++) {
            if(Math.abs(A[i][column]) > Math.abs(A[max][column])) max = i;
        }
        return max;
    }
    private static void swapRows(double[][] A, int i, int k){
        double[] tmp = A[i];
        A[i] = A[k];
        A[k] = tmp;
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

    public static void main(String[] args) {
        double[][] test = new double[][]{{1.0,2.0} , {-2.0,-4.0 }};
        double[] x = solveSing(test);
        System.out.println(Arrays.toString(x));
    }
}
