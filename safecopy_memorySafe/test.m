% %% continue 
% 
% for n = 1:50
%     if mod(n,7)
% %         disp(n);
%         continue
%     end
%     disp(['Divisible by 7: ' num2str(n)])
% end
% 
% 
% %% break
% limit = 0.8;
% s = 0;
% 
% % this jsut lets it run forever until the condition is met and break takes 
% % it out. For could also be used, but then it could happen that the maximum
% % number of iteration is not enough at some point.
% while 1 
%     tmp = rand;
%     if tmp > limit
%         s = s + 1;
%         formatSpec = '%4.2f iterations had to be run with rand command to reach 0.8 value\n';
%         fprintf(formatSpec,s);
%         break
%     end
%     s = s + 1;
% end
% 
% 
%% switch-case-otherwise
% Conclusion1: this is really close to if-elseif-else statements, but
% cleaner view of the code.
% Conclusion2: the best way to make the AND, OR options with the
% switch-case-otherwise setup (at least with my current knowledge) is to
% write down and if-elseif-end statement before and, deppending on the
% conditions, and outcome number (kind of case) is assigned, ie x =
% 1,2,3... that later on gets selected for in the case statements.
% 
% % First option
% n = input('Enter a number: ');
% 
% switch n
%     case -1
%         disp('negative one')
%     case 0
%         disp('zero')
%     case 1
%         disp('positive one')
%     otherwise
%         disp('other value')
% end
% 
% % Second option
% x = NaN;
% switch mod(x, 4)
%     case {0, 3}
%         disp('first')
%     case {1, 2}
%         disp('second')
%     otherwise
%         disp('non-integer')
% end
% 
% % Third option (not really working what i tried)
% x = 1;
% switch x
%     case {1, 2}
%         disp('result is 1 or 2')
%     case {3, 4}
%         disp('result is 1 or 2')
%     case {5, 6}
%         disp('result is 1 or 2')
%     otherwise
%         disp('Result is not 1, 2, 3, 4, 5 or 6.')
% end
% 
% 
%% Try-catch
% try-catch statements run to see the a possible error and further analyze.

% % % Situation 1: It cam be used to give nmre information on the error.
% % The following matrices A and B cannot concatenate
% A = rand(3);
% B = ones(5);
% % C = [A; B];
% 
% % analysis of the error
% try
%    C = [A; B];
% catch ME
%    if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
%       msg = ['Dimension mismatch occurred: First argument has ', ...
%             num2str(size(A,2)),' columns while second has ', ...
%             num2str(size(B,2)),' columns.'];
%         causeException = MException('MATLAB:myCode:dimensions',msg);
%         ME = addCause(ME,causeException);
%    end
%    rethrow(ME)
% end ;
% 
% 
% % % Situation 2: The error can be packaged as a warning and the code then
% % continues.
% try
%     a = notaFunction(5,6);
% catch
%     warning('Problem using function.  Assigning a value of 0.');
%     a = 0;
% end
% disp('The code stil runs');
% % u = 1
% 
% 
% % % Situation 3: Diffrent types of errors can be handled.
% try
%     a = notaFunction(5,6);
% catch ME
%     switch ME.identifier
%         case 'MATLAB:UndefinedFunction'
%             warning('Function is undefined.  Assigning a value of NaN.');
%             a = NaN;
%         case 'MATLAB:scriptNotAFunction'
%             warning(['Attempting to execute script as function. '...
%                 'Running script and assigning output a value of 0.']);
%             notaFunction;
%             a = 0;
%         otherwise
%             rethrow(ME)
%     end
% end
% disp('The code stil runs');
% 
% 
% %%
% figure
% % Construct a figure with subplots and data
% subplot(2,2,1);
% line1 = plot(1:10,rand(1,10),'b');
% title('Axes 1');
% subplot(2,2,2);
% line2 = plot(1:10,rand(1,10),'g');
% title('Axes 2');
% subplot(2,2,3);
% line3 = plot(1:10,rand(1,10),'r');
% title('Axes 3');
% subplot(2,2,4);
% line4 = plot(1:10,rand(1,10),'y');
% title('Axes 4');
% % Construct a Legend with the data from the sub-plots
% hL = legend([line1,line2,line3,line4],{'Data Axes 1','Data Axes 2','Data Axes 3','Data Axes 4'});
% % Programatically move the Legend
% newPosition = [0.4 0.4 0.2 0.2];
% newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits);
% 
% 

%%
fun = @testCase;
% x0 = 0.000000001;
% x0 = 0.9999;
x0 = 0.000149989;

x0array = [1E-6 1E-5 1E-4 1E-3 1E-2 1E-1 1E0 1E1 1E2 1E3];
xres = zeros(1,10);
fvalres = zeros(1,10);
for i = 1:10
%     disp(i);
    x0 = x0array(i);
    [x,fval,exitflag,output] = fsolve(fun,x0);
    xres(i) = x;
    fvalres(i) = fval;
end
% goos solutions
xres(1:3) = 1.4989e-04;
xres(4) = 0.0010;
xres(5) = 0.0042;
% disp(x);

function F = testCase(x)
F = 0.00149*1724 - ((0.15e-3 - x ) / x ) * 1/(0.050 + x) * ((5e-3 - x) / x )*((1e-3 - x ) / x );
% F(1) = 0.00149*1724 - ((0.15e-3 - x(1) ) / x(1) ) * 1/(0.050 + x(1)) * ((5e-3 - x(1)) / x(1) )*((1e-3 - x(1) ) / x(1) );
% % F(1) = 0.00149*1724 - ((0.15e-3 - x(1) ) / x(1) ) * ((5e-3 - x(1)) / x(1) )*((1e-3 - x(1) ) / x(1) );
% % % % F(1) = x / (1-x) * x / (1-x);
% F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
% F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;
end




% fun = @root2d;
% x0 = [0,0];
% x = fsolve(fun,x0);
% 
% function F = root2d(x)
% F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
% F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;
% end