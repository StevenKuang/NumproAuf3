import java.util.Arrays;
import java.util.Comparator;

public class PageRank {

    /**
     * Diese Methode erstellt die Matrix A~ fuer das PageRank-Verfahren
     * PARAMETER:
     * L: die Linkmatrix (s. Aufgabenblatt)
     * rho: Wahrscheinlichkeit, anstatt einem Link zu folgen,
     *      zufaellig irgendeine Seite zu besuchen
     */
    // public static double[][] buildProbabilityMatrix(int[][] L, double rho) {
    //     int length = L.length;
    //     double [][] result = new double[length][length];
    //     //wrong it doesn't work with streams :D we need loops
    //     //Arrays.stream(L).forEach(x-> Arrays.stream(x).map(p-> (int)((1-rho)*(p) + rho/length)).toArray());

    //     for(int i = 0; i < length ; i++) {
    //         for (int j = 0;j < length; j++) {
    //             result[i][j]= ((1 - rho) * (L[i][j]) + rho/length);
    //         }
    //     }
    //     return result;
    // }

    // cuz your version and mine are literally the same thing so i just merged them xD.
    // also we don't need that much loops this way :)
    public static double[][] buildProbabilityMatrix(int[][] L, double rho) {
        
        int dim = L.length;
        double[][] prob = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) 
                sum += L[j][i];
            for (int j = 0; j < dim; j++) {
                if(L[j][i] == 1)
                    prob[j][i] = 1 / sum;
                else
                    prob[j][i] = 0.0;
                prob[j][i] = (1 - rho) * prob[j][i] + rho/dim;
                // 1/4 + 1/6 = 5/12
                // 1/6
                // 5/12
            }
        }
        return prob;
    }


    /**
     * Diese Methode berechnet die PageRanks der einzelnen Seiten,
     * also das Gleichgewicht der Aufenthaltswahrscheinlichkeiten.
     * (Entspricht dem p-Strich aus der Angabe)
     * Die Ausgabe muss dazu noch normiert sein.
     * PARAMETER:
     * L: die Linkmatrix (s. Aufgabenblatt)
     * rho: Wahrscheinlichkeit, zufaellig irgendeine Seite zu besuchen
     * ,anstatt einem Link zu folgen.
     *
     */
    public static double[] rank(int[][] L, double rho) {
       double[][] probMatrix = buildProbabilityMatrix(L, rho);
       int length = probMatrix.length;
        for (int i = 0; i < length; i++) {
            probMatrix[i][i] -= 1.0;
        }
        double[] p = Gauss.solveSing(probMatrix);
        double sum = Arrays.stream(p).sum();
        if(sum == 0) return p;
        double one = 0;
        for (int i = 0; i < p.length; i++) {
        p[i] *= (1/sum);
        }
        return p;
    }

    /**
     * Diese Methode erstellt eine Rangliste der uebergebenen URLs nach
     * absteigendem PageRank.
     * PARAMETER:
     * urls: Die URLs der betrachteten Seiten
     * L: die Linkmatrix (s. Aufgabenblatt)
     * rho: Wahrscheinlichkeit, anstatt einem Link zu folgen,
     *      zufaellig irgendeine Seite zu besuchen
     */
    public static String[] getSortedURLs(String[] urls, int[][] L, double rho) {
        int n = L.length;

        double[] p = rank(L, rho);

        RankPair[] sortedPairs = new RankPair[n];
        for (int i = 0; i < n; i++) {
            sortedPairs[i] = new RankPair(urls[i], p[i]);
        }

        Arrays.sort(sortedPairs, new Comparator<RankPair>() {

            @Override
            public int compare(RankPair o1, RankPair o2) {
                return -Double.compare(o1.pr, o2.pr);
            }
        });

        String[] sortedUrls = new String[n];
        for (int i = 0; i < n; i++) {
            sortedUrls[i] = sortedPairs[i].url;
        }

        return sortedUrls;
    }

    /**
     * Ein RankPair besteht aus einer URL und dem zugehoerigen Rang, und dient
     * als Hilfsklasse zum Sortieren der Urls
     */
    private static class RankPair {
        public String url;
        public double pr;

        public RankPair(String u, double p) {
            url = u;
            pr = p;
        }
    }
    // public static void main(String[] args) {
    //     int[][] L = {{1,0,1},{0,1,0},{1,1,1}};
    //     double rho = 0.5;
    //     buildProbabilityMatrix(L, rho);
    //     System.out.println("res");
    // }
}


