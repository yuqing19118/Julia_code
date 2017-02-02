
# Part e goes in this Julia cell. 

# Question 4 code goes in this Julia cell 

Pkg.add("Dierckx")
Pkg.update()
using Dierckx
using PyPlot
############################## Part a ##############################
function routine(n,a,b,x,f)
    #inputs1 holds n equidistant points on the interval [a,b], outputs1 holds corresponding y values
    inputs1 = zeros(n);
    outputs1 = zeros(n);
    #inputs2 holds n Chebyshev points on the interval [a,b]
    inputs2 = zeros(n);
    #interpolant1 is the polynomial interpolant to f at n equidistant points on [a,b] using Newton polynomials
    interpolant1 = 0;
    #interpolant2 is the polynomail interpolant to f at n Chebyshev points on [a,b] using Newton polynomials
    interpolant2 = 0;
    #interpolant3 is the cubic spline interpolant at n equidistant points on [a,b]
    interpolant3 = 0;
    
    for i=1:n
        #Calculate n equidistant points on [a,b] and store to inputs1
        inputs1[i] = a + (b-a)/(n-1)*(i-1)
        outputs1[i] = f(inputs1[i])
        #Calculate n Chebyshev points on [a,b] and store to inputs2
        inputs2[i] = (a+b)/2 + (b-a)/2*cos((i-1)*pi/(n-1))
    end
    
    #divided1 holds divided difference using n equidistant points on [a,b]
    divided1 = zeros(n,n)
    #divided2 holds divided difference using n Chebyshev points on [a,b]
    divided2 = zeros(n,n)
   
    #start from first column and last row, calculate the divided difference
    for column = 1:n
        for row = n+1-column:-1:1
            if column == 1
                #if we are at the first column, then the divided difference is f(x[row])
                divided1[row,column] = f(inputs1[row])    
                divided2[row,column] = f(inputs2[row])
            else
                #else calculate according to the formula given in lecture notes
                divided1[row,column] = (divided1[row,column-1]-divided1[row+1,column-1])/(inputs1[row]-inputs1[row+column-1])
                divided2[row,column] = (divided2[row,column-1]-divided2[row+1,column-1])/(inputs2[row]-inputs2[row+column-1])
            end
        end
   end
   
    #Print the divided difference of a polynomial interpolant at n equidistant points on [a,b]
    println("divided difference of a polynomial interpolant at ",n," equidistant points on [a,b] ********************")
    println(divided1)
    
    #Print the divided difference of a polynomial interpolant at n Chebyshev points on [a,b]
    println("divided difference of a polynomial interpolant at ",n," Chebyshev points on [a,b] *********************")
    println(divided2)

    #the first term in the interpolants of 1st and 2nd interpolation is D[0,0]*1
    interpolant1 = divided1[1,1];
    interpolant2 = divided2[1,1];
    #start from the second term in the interpolant
    for j=1:n-1
        newton1 = 1
        newton2 = 1
        for k=1:j
            newton1 =  newton1 .* (x-inputs1[k])
            newton2 = newton2.*(x-inputs2[k])
        end
        interpolant1 = interpolant1 + divided1[1,j+1].*newton1
        interpolant2 = interpolant2 + divided2[1,j+1].*newton2
    end
    
    #Calculate cubic spline interpolant at n equidistant points on [a,b] using the sample code given
    spl = Spline1D(inputs1, outputs1)
    interpolant3 = evaluate(spl,x)
    
    return interpolant1, interpolant2, interpolant3
end



############################## Part b ##############################
f1 = x -> exp(x)
f2 = x -> sin(x)
f3 = x -> log(x+3)
f4 = x -> abs(x).^(3/2)
a = 0
b = 1


println("---------- Outputs for the first function f1 = exp(x) at n = 5 ------------------------------------------------")
n = 5
x = linspace(a,b,n)
result1, result2, result3 = routine(n,a,b,x,f1)
title("Graph of error on the interval [0,1] using interpolants in 4a with n = 5 on the function f1 = exp(x)")
plot(x, abs(f1(x)-result1),color = "red", label = "polynomial interpolant equidistant points error")
plot(x, abs(f1(x)-result2), color = "green", label = "polynomial interpolant Chebyshev points error")
plot(x, abs(f1(x)-result3), color = "blue",label = "cubic spline error")
legend()
println("****Observation****")
println("The polynomial interpolation using Chebyshev points is worse than the polynomial", 
"using equidistant points and cubic spline")


println("---------- Outputs for the first function f1 = exp(x) at n = 12 ------------------------------------------------")
n = 12
x = linspace(a,b,n)
result1, result2, result3 = routine(n,a,b,x,f1)
title("Graph of error on the interval [0,1] using interpolants in 4a with n = 12 on the function f1 = exp(x)")
plot(x, abs(f1(x)-result1),color = "red", label = "polynomial interpolant equidistant points error")
plot(x, abs(f1(x)-result2), color = "green", label = "polynomial interpolant Chebyshev points error")
plot(x, abs(f1(x)-result3), color = "blue",label = "cubic spline error")
legend()
println("****Observation****")
println("The cubic spline interpolant is better on a whole than each of the polynomial interpolants. While",
" polynomial interpolant using equidistant points performs well on the first half of the interval, the average error ",
"using cubic spline seems smaller than each of the polynomial interpolants")

println("---------- Outputs for the first function f2 = sin(x) at n = 5 ------------------------------------------------")
n = 5
x = linspace(a,b,n)
result1, result2, result3 = routine(n,a,b,x,f2)
title("Graph of error on the interval [0,1] using interpolants in 4a with n = 5 on the function f2 = sin(x)")
plot(x, abs(f2(x)-result1),color = "red", label = "polynomial interpolant equidistant points error")
plot(x, abs(f2(x)-result2), color = "green", label = "polynomial interpolant Chebyshev points error")
plot(x, abs(f2(x)-result3), color = "blue",label = "cubic spline error")
legend()
println("****Observation****")
println("The polynomial interpolation using Chebyshev points is a little bit worse than the polynomial", 
"using equidistant points and cubic spline")
         

println("---------- Outputs for the first function f2 = sin(x) at n = 12 ------------------------------------------------")
n = 12
x = linspace(a,b,n)
result1, result2, result3 = routine(n,a,b,x,f2)
title("Graph of error on the interval [0,1] using interpolants in 4a with n = 12 on the function f2 = sin(x)")
plot(x, abs(f2(x)-result1),color = "red", label = "polynomial interpolant equidistant points error")
plot(x, abs(f2(x)-result2), color = "green", label = "polynomial interpolant Chebyshev points error")
plot(x, abs(f2(x)-result3), color = "blue",label = "cubic spline error")
legend()
println("****Observation****")
println("The polynomial interpolation using equidistant points is better than the other two on first half of the interval. ", 
"On the second half of the interval, cubic spline seems better, but all three interpolants are not accurate. Overall, ",
"cubic spline interpolant has lower average error.")

println("---------- Outputs for the first function f3 = log(x+3) at n = 5 ------------------------------------------------")
n = 5
x = linspace(a,b,n)
result1, result2, result3 = routine(n,a,b,x,f3)
title("Graph of error on the interval [0,1] using interpolants in 4a with n = 5 on the function f3 = log(x+3)")
plot(x, abs(f3(x)-result1),color = "red", label = "polynomial interpolant equidistant points error")
plot(x, abs(f3(x)-result2), color = "green", label = "polynomial interpolant Chebyshev points error")
plot(x, abs(f3(x)-result3), color = "blue",label = "cubic spline error")
legend()
println("****Observation****")
println("The polynomial interpolation using Chebyshev points is worse than the polynomial", 
"using equidistant points and cubic spline")

println("---------- Outputs for the first function f3 = log(x+3) at n = 15 ------------------------------------------------")
n = 15
x = linspace(a,b,n)
result1, result2, result3 = routine(n,a,b,x,f3)
title("Graph of error on the interval [0,1] using interpolants in 4a with n = 15 on the function f3 = log(x+3)")
plot(x, abs(f3(x)-result1),color = "red", label = "polynomial interpolant equidistant points error")
plot(x, abs(f3(x)-result2), color = "green", label = "polynomial interpolant Chebyshev points error")
plot(x, abs(f3(x)-result3), color = "blue",label = "cubic spline error")
legend()
println("****Observation****")
println("The polynomial interpolation using equidistant points is better than the other two interpolants. ", 
"What is notable is that the error of polynomial interpolation using Chebyshev points decreases as x increases.")

println("---------- Outputs for the first function f4 = abs(x)^(3/2) at n = 6 ------------------------------------------------")
n = 6
x = linspace(a,b,n)
result1, result2, result3 = routine(n,a,b,x,f4)
title("Graph of error on the interval [0,1] using interpolants in 4a with n = 6 on the function f4 = abs(x)^(3/2)")
plot(x, abs(f4(x)-result1),color = "red", label = "polynomial interpolant equidistant points error")
plot(x, abs(f4(x)-result2), color = "green", label = "polynomial interpolant Chebyshev points error")
plot(x, abs(f4(x)-result3), color = "blue",label = "cubic spline error")
legend()
println("****Observation****")
println("The polynomial interpolation using Chebyshev points is worse than the polynomial", 
"using equidistant points and cubic spline", " but the error of polynomial interpolant using Chebyshev points increases first ",
"then decreases as x increases.")

println("---------- Outputs for the first function f4 = abs(x)^(3/2) at n = 12 ------------------------------------------------")
n = 12
x = linspace(a,b,n)
result1, result2, result3 = routine(n,a,b,x,f4)
title("Graph of error on the interval [0,1] using interpolants in 4a with n = 12 on the function f4 = abs(x)^(3/2)")
plot(x, abs(f4(x)-result1),color = "red", label = "polynomial interpolant equidistant points error")
plot(x, abs(f4(x)-result2), color = "green", label = "polynomial interpolant Chebyshev points error")
plot(x, abs(f4(x)-result3), color = "blue",label = "cubic spline error")
legend()
println("****Observation****")
println("The polynomial interpolation using Chebyshev points is a little worse than the polynomial", 
"using equidistant points and cubic spline", " and the error of polynomial interpolant using Chebyshev points oscillates ",
"but displays a decreasing trend. The error using 12 Chebyshev points ",
"is significantly lower than using 6 Chebyshev points.")


