package difequationsapp;

import static java.lang.Math.exp;


public class DifEquationsApp {


    public static void main(String[] args) {
        int a = 0;
        int b = 1;
        int n = 10;
        double eps = 0.000001;
        double h = ((double)b - a) / n;
        int size = n + 1;
        
        eulerMod(a, b, n, h, size, eps);
        rungeCutta2(a, b, n, h, size);
        rungeCutta4(a, b, n, h, size);
        adams(a, b, n, h, size);

    }
    
    static void rungeCutta2(int a, int b, int n, double h0, int size){
        double[] x = new double[size];
        x[0] = (double)a;
        for(int i = 1; i < size; i++)
            x[i] = x[0] + i*h0;
        
        double k1, k2;
        
        double[] y = new double[size];
        y[0] = 0;
       
        for(int i = 1; i < size; i++){
            k1 = h0 * funcValue(x[i-1], y[i-1]);
            k2 = h0 * funcValue(x[i-1] + h0, k1);
            y[i] = y[i-1] + (k1 + k2) / 2.0; 
        }
        System.out.println("");
        System.out.println("y' = xe^(-x^2) - 2xy");
        System.out.println("y(0) = 0");
        System.out.println("[" + a + "][" + b +"]");
        System.out.println("Метод Рунге-Кутты(2 порядка точности)");
        
        System.out.println("Решение: ");
        System.out.print("x: ");
        printArr(x, size);
        System.out.print("y: ");       
        printArr(y, size);
        System.out.print("y(т): ");
        for(int i = 0; i < size; i++)
            System.out.print(solutionValue(x[i]) + " ");
        System.out.println("\nПогрешности:");
        for(int i = 0; i < size; i++)
            System.out.print(Math.abs(solutionValue(x[i]) - y[i]) + " ");
        System.out.println("");
    }
    
    static void rungeCutta4(int a, int b, int n, double h0, int size){
        double[] x = new double[size];
        x[0] = (double)a;
        for(int i = 1; i < size; i++)
            x[i] = x[0] + i*h0;
        
        double k1, k2, k3, k4;
        
        double[] y = new double[size];
        y[0] = 0;
       
        for(int i = 1; i < size; i++){
            k1 = h0 * funcValue(x[i-1], y[i-1]);
            k2 = h0 * funcValue(x[i-1]+h0 / 2.0, y[i-1] + k1 / 2.0);
            k3 = h0 * funcValue(x[i-1]+h0 / 2.0, y[i-1] + k2 / 2.0);
            k4 = h0 * funcValue(x[i-1]+h0, y[i-1] + k3);
            y[i] = y[i-1] + (k1 + 2*(k2 + k3) + k4) / 6.0; 
        }
        System.out.println("");
        System.out.println("y' = xe^(-x^2) - 2xy");
        System.out.println("y(0) = 0");
        System.out.println("[" + a + "][" + b +"]");
        System.out.println("Метод Рунге-Кутты(4 порядка точности)");
        
        System.out.println("Решение: ");
        System.out.print("x: ");
        printArr(x, size);
        System.out.print("y: ");       
        printArr(y, size);
        System.out.print("y(т): ");
        for(int i = 0; i < size; i++)
            System.out.print(solutionValue(x[i]) + " ");
        System.out.println("\nПогрешности:");
        for(int i = 0; i < size; i++)
            System.out.print(Math.abs(solutionValue(x[i]) - y[i]) + " ");
        System.out.println("");
    }
    
    static void adams(int a, int b, int n, double h, int size){
        double[] x = new double[size];
        x[0] = a;
        for(int i = 1; i < size; i++)
            x[i] = x[0] + i*h;
        
        double k1, k2, k3, k4;
        
        double[] y = new double[size];
        y[0] = 0;
       
        for(int i = 1; i < 4; i++){
            k1 = h * funcValue(x[i-1], y[i-1]);
            k2 = h * funcValue(x[i-1]+h / 2.0, y[i-1] + k1 / 2.0);
            k3 = h * funcValue(x[i-1]+h / 2.0, y[i-1] + k2 / 2.0);
            k4 = h * funcValue(x[i-1]+h, y[i-1] + k3);
            y[i] = y[i-1] + (k1 + 2*(k2 + k3) + k4) / 6.0; 
        }
        
        double[] yP = new double[size];
        double[] fP = new double[size];
        for(int i = 3; i < size-1; i++){
            yP[i+1] = y[i] + h/24.0 * (55*funcValue(x[i], y[i]) - 59*funcValue(x[i-1], y[i-1]) + 37*funcValue(x[i-2], y[i-2]) - 9*funcValue(x[i-3], y[i-3]));
            fP[i+1] = funcValue(x[i+1], yP[i+1]);
            y[i+1] = y[i] + h/24.0 * (9*fP[i+1] + 19*funcValue(x[i], y[i]) - 5*funcValue(x[i-1], y[i-1]) + funcValue(x[i-2], y[i-2]));
        }
        
        System.out.println("");
        System.out.println("y' = xe^(-x^2) - 2xy");
        System.out.println("y(0) = 0");
        System.out.println("[" + a + "][" + b +"]");
        System.out.println("Метод Адамса");
        System.out.println("Решение: ");
        System.out.print("x: ");
        printArr(x, size);
        System.out.print("y: ");       
        printArr(y, size);
        System.out.print("y(т): ");
        for(int i = 0; i < size; i++)
            System.out.print(solutionValue(x[i]) + " ");
        System.out.println("\nПогрешности:");
        for(int i = 0; i < size; i++)
            System.out.print(Math.abs(solutionValue(x[i]) - y[i]) + " ");
        System.out.println("");
    }
    
    static void eulerMod(int a, int b, int n, double h, int size, double eps){
        double[] x = new double[size];
        x[0] = (double)a;
        for(int i = 1; i < size; i++)
            x[i] = x[0] + i*h;
        
        double[] y = new double[size];
        y[0] = 0;
        for(int i = 1; i < size; i++)
            y[i] = y[i-1] + h*funcValue(x[i-1], y[i-1]);
        
        System.out.println("y' = xe^(-x^2) - 2xy");
        System.out.println("y(0) = 0");
        System.out.println("[" + a + "][" + b +"]");
        System.out.println("Метод Эйлера(модифицированный)");
        
        System.out.println("Решение: ");
        System.out.print("x: ");
        printArr(x, size);
        System.out.print("y: ");
        printArr(y, size);
        
        while(Math.abs(y[1] - solutionValue(x[1])) >= eps)                      //Уточнение решения
            for(int i = 1; i < size; i++)
                y[i] = y[i-1] + h/2 * (funcValue(x[i-1], y[i-1]) + funcValue(x[i], y[i]));
           
        System.out.print("y*: ");
        printArr(y, size);
        
        System.out.print("y(т): ");
        for(int i = 0; i < size; i++)
            System.out.print(solutionValue(x[i]) + " ");
        System.out.println("\nПогрешности:");
        for(int i = 0; i < size; i++)
            System.out.print(Math.abs(solutionValue(x[i]) - y[i]) + " ");
        System.out.println("");
    }
    
    static double funcValue(double x, double y){
        return x * exp(-Math.pow(x, 2)) - 2 * x * y;
    }
    
    static double solutionValue(double x){
        return exp(-Math.pow(x, 2)) * Math.pow(x, 2) / 2.0;
    }
    
    static void printArr(double[] arr, int size){
        for(int i = 0; i < size; i++)
            System.out.print(arr[i] + "\t");
        System.out.println();
    }
    
}
