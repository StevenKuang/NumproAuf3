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
        return backSubst(A, b);
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
        double[][] copyA = new double[A.length][A.length];
        for (int i = 0; i < A.length; i++) {
            copyA[i] = Arrays.copyOf(A[i], A[i].length);
        }
        int pivot = 0; double coefficient = 1.0; boolean noPivot= false; int k;
        for (k = 0; k < A.length; k++) {
            pivot = getIndexRowPivotUnderElement(copyA,k,k);
            if(A[pivot][k] < 1E-10) { noPivot = true; break;}
            swapRows(copyA,k,pivot);
            // Gauss
            for (int i = k + 1; i < A.length; i++) {
                coefficient = copyA[i][k] / copyA[k][k];
                for (int j = k; j < A.length; j++) {
                    copyA[i][j] -= coefficient * copyA[k][j];
                }
            }
        }
        if(noPivot) {
            // extract the vector v which the column where we stopped the gauss elimination => k
            double[] v = new double[k+1];
            for (int i = 0; i <k+1  ; i++) {
                v[i] = - copyA[i][k];
            }
            // extract the matrix which is the part of the original matrix from 0->k-1
            double[][] T = new double[k+1][k+1];
            for (int i = 0; i < k+1; i++) {
                System.arraycopy(copyA[i], 0, T[i], 0, k+1);
            }
            double[] x = new double[A.length];
            double sum;
            for (int i = T.length - 1; i >= 0; i--) {
                sum = 0.0;
                for (int j = i + 1; j < T.length; j++) {
                    sum += T[i][j] * x[j];
                }
                x[i] = (v[i] - sum) / T[i][i];
            }
            // the last index that was filled was k
            x[k+1] = 1.0;
            for (int i = k+2; i < x.length; i++) {
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
