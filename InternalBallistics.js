
// case capacity source:  http://kwk.us/cases.html
// 253 grains of water per cubic inch 

// various global state flags
var bestSet;
var keepRunning;
var generations;
var targetError;
var loadIndex;

// target values for a given simulation iteration
var targetPSI;
var targetFPS;

// this will be an array of arrays; first index is the individual's number,
//  the second index will be the parameters. 
var parameterSet;

// set of parameters for powder definition
const parmCount = 4;
var pRange = [];
var pMin = [];

var grainsToSlugs = 1 / 7000 / 32.2;

// IMR 3031 55gr 24"
// 20.0gr 2878fps 46.8kpsi
var charge = 20.0;
var bulletMass = 55 * grainsToSlugs;
var bore = 0.223;
var boreArea = Math.pow(bore / 2, 2) * 3.14;
var barrelLen = 24;
var chamberVolume = 0.112;    // initial volume of case, 1.5" long times double bore area
var dragStatic = 250;    // resistance to movements in pounds
var dragKinetic = 100;

// for graphing pressure and velocity
var xscale; // fps max
var yscale;  // kpsi max

function Init()
{
    var outputDoc = document.getElementById("outputTextArea");
    outputDoc.innerHTML += "In Init function.\n";

	loadIndex = 0;

    // array of 100 sest of 3 parameters, plus the error
    parameterSet = [];
    for(var i = 0; i < 100; ++i)
    {
        parameterSet[i] = [];
    }

    bestSet = [];
    bestSet[0] = 0;
    bestSet[1] = 0;
    bestSet[2] = 0;
    bestSet[3] = 0;
    bestSet[parmCount] = -1;  // negative error means it hasn't been set yet

    for(var set = 0; set < 100; ++set)
    {
        parameterSet[set][0] = pMin[0] + Math.random() * pRange[0];
        parameterSet[set][1] = pMin[1] + Math.random() * pRange[1];
        parameterSet[set][2] = pMin[2] + Math.random() * pRange[2];
        parameterSet[set][3] = pMin[3] + Math.random() * pRange[3];

        parameterSet[set][parmCount] = -1;
    }

	// minimum and range of each parameter in powder definition
	pRange[0] = 1;  // curve shape parameter, 0-1
	pMin[0] = 0;
	pRange[1] = 10;  // curve slope parameter
	pMin[1] = 0;
	pRange[2] = 10;  // pressure modifier
	pMin[2] = 0;
	pRange[3] = 2000;  // gas volume per grain of powder, cubic inches per grain
	pMin[3] = 100;

	xscale = 750 / 5000; // fps max
	yscale = 500 / 100;  // kpsi max
}

// called when error is less than target, or user clicks stop
function StopGenetic()
{
    var outputDoc = document.getElementById("outputTextArea");
    keepRunning = false;

    outputDoc.innerHTML = "Best error: " + bestSet[parmCount] + "\n";
    outputDoc.innerHTML += "Parameter values:\n";
    // generate 90 new individuals by mixing and mutating top 10
    for(var i = 0; i < parmCount; ++i)
    {
        outputDoc.innerHTML += "Parameter " + i + ":  " + bestSet[i] + "\n";
    }
   
    outputDoc.innerHTML += "Generations needed, " + generations + "\n";
}

// simulation iteration loop
function Iterate()
{
    var outputDoc = document.getElementById("outputTextArea");
    timeStep = 10; // ms

    var hitTarget = 0;
    targetError = 1;
	currentError = 100;

    if(keepRunning)
    {
		// set up load data
		charge = loads[loadIndex][0];
		bulletMass = loads[loadIndex][1] * grainsToSlugs;
		bore = loads[loadIndex][2];
		boreArea = Math.pow(bore / 2, 2) * Math.pi;
		barrelLen = loads[loadIndex][3];
		chamberVolume = loads[loadIndex][4];
		dragStatic = loads[loadIndex][5];
		dragKinetic = loads[loadIndex][6];
	    targetFPS = loads[loadIndex][7];
    	targetPSI = loads[loadIndex][8];

		loadIndex += 1;

		// wrap around when we hit the end of the list
		if(loadIndex >= loads.length) loadIndex = 0;

        Draw();
        // if we have a best set, and it's really close, terminate
        if(bestSet[parmCount] > 0 && bestSet[parmCount] < currentError)
        {
			// work acceptable error down to target slowly, so we
			//  don't over-optimize on one solution at first
			if(currentError < targetError) --currentError;

			// if we do hit the target, increment a counter and keep
			//  going until we hit it repeatedly
			if(bestSet[parmCount] < targetError)
				++hitTarget;

			if(hitTarget > 25)
	            StopGenetic();

            return;
        }
        Evolve();
   
        var timerID = setTimeout(Iterate, timeStep);
    }
}

// copy parameters in set a to set b
function Clone(a, b)
{
    for(var i = 0; i < parmCount + 1; ++i)
    {
        parameterSet[b][i] = parameterSet[a][i];
    }
}

// copy parameters from best into set a
function CloneBest(a)
{
    for(var i = 0; i < parmCount; ++i)
    {
        parameterSet[a][i] = bestSet[i];
    }
}

// randomly swap parameters between sets a and b
function Breed(a, b)
{
    var outputDoc = document.getElementById("outputTextArea");

    // remember old individuals, so we can ensure new ones coming out
    var oldA = a;
    var oldB = b;

    // keep swapping genes until we get two new individuals, or give up
    var attempt = 0;
    var success = false;
    while(attempt < 10)
    {
        ++attempt;
        // go through each parameter
        for(var i = 0; i < parmCount; ++i)
        {
            // make sure the genes we're flipping differ
            if(parameterSet[a][i] != parameterSet[b][i])
            {
                // flip a coin to decide if these parameters swap
                if(Math.random() > 0.5)
                {
                    var temp = parameterSet[a][i];
                    parameterSet[a][i] = parameterSet[b][i];
                    parameterSet[b][i] = temp;
                }
            }
        }
        if(a != oldA && a != oldB && b != oldA && b != oldB)
        {
            success = true;
            break;
        }
    }

    // if breeding was unsuccessful, mutate instead
    if(!success)
    {
        outputDoc.innerHTML += "Breeding failed, mutating instead.\n";
        Mutate(a);
        Mutate(b);
    }
}

// randomly adjust parameters by a smallish amount
function Mutate(a)
{
    // go through each parameter
    for(var i = 0; i < parmCount; ++i)
    {
        // flip a coin to decide if the parameter changes
        var flip = Math.random();
        // 50% chance of a small mutation, +/- 10%
        if(flip > 0.5)
        {
            parameterSet[a][i] += ((Math.random() - 0.5) * pRange[i] / 10);
        }
        // 10% chance of a big mutation, +/- 50%
        if(flip < 0.1)
        {
            parameterSet[a][i] += ((Math.random() - 0.5) * pRange[i] / 2);
        }
        // clip mutation to parameter range
        if(parameterSet[a][i] > pMin + pRange[i]) parameterSet[a][i] = pMin + pRange[i];
        if(parameterSet[a][i] < pMin) parameterSet[a][i] = pMin;
    }
}

// function to sort lowest error into index 0 of array
//  error is last parameter in set
function compareError(a, b)
{
    // return less than 0 if a < b
    return a[parmCount] - b[parmCount];
}

// keep the fittest, kill the rest, breed and mutate the survivors to make more
function Evolve()
{
    var outputDoc = document.getElementById("outputTextArea");
    ++generations;
    // sort the individuals by error, put the best in the low indicies
    parameterSet.sort(compareError);

    // see if we've beat the best set; if we haven't picked one, index 0 is best
    if(bestSet[parmCount] < 0 || parameterSet[0][parmCount] < bestSet[parmCount])
    {
        for(var i = 0; i < parmCount + 1; ++i)
        {
            bestSet[i] = parameterSet[0][i];
        }
    }

	outputDoc.innerHTML = "Target " + targetPSI + "kpsi and " + targetFPS + "fps\n";

    outputDoc.innerHTML += "Best error: " + bestSet[parmCount] + "\n";
    outputDoc.innerHTML += "Generation " + generations + " survivors:\n";
    // generate 90 new individuals by mixing and mutating top 10
    outputDoc.innerHTML += "Error: " + parameterSet[0][parmCount] + "\n";
    for(var i = 1; i < 10; ++i)
    {
        outputDoc.innerHTML += "Error: " + parameterSet[i][parmCount];
        if(
            parameterSet[i][3] != parameterSet[0][3] ||
            parameterSet[i][2] != parameterSet[0][2] ||
            parameterSet[i][1] != parameterSet[0][1] ||
            parameterSet[i][0] != parameterSet[0][0]
        )
        {
            outputDoc.innerHTML += "\n";
        }
        else
        {
            outputDoc.innerHTML += ", clone\n";
        }
    }

    // start from bottom, fill up to 0 so we keep the top 10
    // 5 mutants of best
    outputDoc.innerHTML += "5 mutants of best set\n";
    var i = 99;
    for(var j = 0; j < 5; ++j)
    {
        CloneBest(i);
        Mutate(i);
        --i;
    }
    // i should be 94
    // 20 breeds of best with top 10
    outputDoc.innerHTML += "20 hybrids of best and each of top 10\n";
    for(var j = 0; j < 10; ++j)
    {
        CloneBest(i);
        Clone(j, i - 1);
        Breed(i, i - 1);
        i -= 2;
    }
    // i should be 74
    // 50 breeds of random top 10
    outputDoc.innerHTML += "50 hybrids between random pairs of top 10\n";
    for(var j = 0; j < 25; ++j)
    {
        // pick 2 different sources
        var source1 = Math.trunc(Math.random() * 10);
        var source2 = source1;
           while(source2 == source1) source2 = Math.trunc(Math.random() * 10);

        Clone(source1, i);
        Clone(source2, i - 1);
        Breed(i, i - 1);
        i -= 2;
    }
    // i should be 24
    // 5 new random individuals
    outputDoc.innerHTML += "5 completely new random individuals\n";
    for(var j = 0; j < 5; ++j)
    {
        parameterSet[i][0] = pMin[0] + Math.random() * pRange[0];
        parameterSet[i][1] = pMin[1] + Math.random() * pRange[1];
        parameterSet[i][2] = pMin[2] + Math.random() * pRange[2];
        parameterSet[i][3] = pMin[3] + Math.random() * pRange[3];
    }
    // i should be 19
    // 20 mutants of top 10
    outputDoc.innerHTML += "20 mutations of top 10 individuals\n";
    for(var j = 0; j < 10; ++j)
    {
        // copy top 10 element into element 10 down
        Clone(i - 10, i);
        // mutate daughter
        Mutate(i);
        --i;
    }
    // i should be 9
    for(var j = 0; j < 10; ++j)
    {
        Mutate(i);
        --i;
    }
    // i should be -1
}

// parameter0 and percent between 0 and 1
//  returns a slope at the given point on the gamma curve
function GammaDelta(percent, parameter0)
{
    // clip input to just inside 0-1 range
    if(percent < 1) percent = 1;
    if(percent > 99) percent = 99;

    // set so that .50 generates a flat line
    parameter0 = parameter0 * 2;
    // cube > 1 values so we get a more symmetric set of curves
    if(parameter0 > 1) parameter1 = Math.pow(parameter0, 3);

    // do basic gamma curve in two close spots
    var val1 = Math.pow(percent - 0.5, parameter0);
    var val2 = Math.pow(percent + 0.5, parameter0);

    // calculate slope
    slope = (val2 - val1);

    return slope;
}

// percent is percent of powder burned, pressure is in kpsi
function BurnPowder(percent, pressure, dt, parameter0, parameter1, parameter2)
{
    // get base burn rate
    var baseRate = GammaDelta(percent, parameter0);

    // adjust slope based on parameter1, with 1 being no change
    baseRate *= parameter1;

    // adjust for pressure, divided by 100 kpsi to keep values in the single digits
    baseRate *= (1 + pressure / 100 * parameter2);

    // fudge factor adjusts burn rate to reasonable order of magnitude
    percent += baseRate * dt * 500;

    if(percent > 1.0) percent = 1.0;
    return percent;
}

// run the simulation using parameter set i, return peak pressure and velocity
function RunSimulation(i, results)
{
    var outputCanvas = document.getElementById("drawArea");
    var gc = outputCanvas.getContext("2d");
    var outputDoc = document.getElementById("outputTextArea");
	var status = document.getElementById("status");

    gc.strokeStyle = "gray";
   
    var percent = 0;    // amount of powder burned
    var velocity = 0;    // inches per second
    var position = 0;    // inches, position of base of bullet
    var volume = chamberVolume;
    var gas = 1 * volume;    // set initial pressure to 1 kpsi
    var pressure = gas / volume / 1000; // kpsi
    var maxPressure = 0;

    var dt = 0.00001;

   	status.innerHTML = "Beginning simulation";

    gc.moveTo(0, 500);
    for(var t = 0; t < 0.1; t += dt)
    {
   		status.innerHTML = "t = " + t;
        // calculate powder burned
        percent = BurnPowder(percent, pressure, dt,
                        parameterSet[i][0],
                        parameterSet[i][1],
                        parameterSet[i][2]);

        // calculate new gas volume
        gas = percent * charge * parameterSet[i][3];

        // calculate pressure in kpsi based on volume and gas amount
        pressure = gas / volume / 1000;

        //outputDoc.innerHTML += "t " + t + " % " + percent + " kpsi " + pressure + "\n";
       
        if(pressure > maxPressure) maxPressure = pressure;

        // bail early if we exceeded 120% of target pressure
        if(pressure > 1.2 * targetPSI)
        {
    		results.psi = pressure;    // kpsi
		    results.fps = velocity / 12;         // inches per second to feet per second

        	status.innerHTML = "Failure: " + result.psi + "kpsi " + result.fps + "fps";
            return;
        }

        // calculate force in pounds on the bullet, remember pressure is kpsi
        var force = pressure * 1000 * boreArea;
        if(velocity > 0)
            force -= dragStatic;
        else
            force -= dragKinetic;

        if(force < 0) force = 0;
        // calculate acceleration of bullet;
        var acceleration = force / bulletMass  * 12; // pounds to inches per second squared
        velocity += acceleration * dt; // inches per second
        position += velocity * dt; // calculate displacement
        volume = position * boreArea + chamberVolume; // increase swept volume

        // bail if we hit the end of the barrel
        if(position > barrelLen) 
		{
    		results.psi = pressure;    // kpsi
		    results.fps = velocity / 12;         // inches per second to feet per second

        	status.innerHTML = "Success: " + result.psi + "kpsi " + result.fps + "fps";
            return;
		}

        // plot the pressure/velocity curve of the best individual
        if(i == 0)
        {
            //gc.lineTo(velocity / 12 * xscale, pressure * yscale);
            //gc.lineTo(position / barrelLen * 750, 500 - pressure * yscale);
            gc.lineTo(t / 0.001 * 750, 500 - pressure * yscale);
            gc.stroke();
        }
    }

    results.psi = maxPressure;    // kpsi
    results.fps = velocity / 12;         // inches per second to feet per second

   	status.innerHTML = "Failure: " + result.psi + "kpsi " + result.fps + "fps";
	return;
}

function RunBallistics()
{
    var outputDoc = document.getElementById("outputTextArea");
    outputDoc.innerHTML = "Running ballistics.\n";

    Init();

    generations = 0;
    keepRunning = true;

    Iterate();
}

function RunTest()
{
    var result = new Object();
    result = {psi:0, fps:0};

    var outputDoc = document.getElementById("outputTextArea");

    parameterSet = [];
    for(var i = 0; i < 100; ++i)
    {
        parameterSet[i] = [];
    }
    parameterSet[0][0] = 0.5;    // curvature
    parameterSet[0][1] = 5;    // slope
    parameterSet[0][2] = 1;    // pressure modifier
    parameterSet[0][3] = 500;    // gas volume in^3 per grain

    outputDoc.innerHTML += "\n\nRunning single test\n";
    RunSimulation(0, result);
    outputDoc.innerHTML += result.psi + "kpsi " + result.fps + "fps\n";
}

function Draw()
{
    var outputDoc = document.getElementById("outputTextArea");
    var outputCanvas = document.getElementById("drawArea");
    var gc = outputCanvas.getContext("2d");
    // clear and draw target
    gc.fillStyle = "white";
    gc.fillRect(0, 0, 750, 500);

    gc.fillStyle = "red";
    gc.fillRect(targetFPS * xscale, 500 - targetPSI * yscale, 10, 10);

    // plot points and calculate error
    for(i = 0; i < 100; ++i)
    {
        var result = new Object();
        result = {psi:0, fps:0};

        RunSimulation(i, result);

        // squared error works, since all we care about is who's closer
        parameterSet[i][parmCount] = Math.pow(targetPSI - result.psi, 2) + Math.pow(targetFPS - result.fps, 2);

        var x = result.fps * xscale;
        var y = 500 - result.psi * yscale;

        // clip x and y to range
        if(x < 1) x = 1;
        if(x > 748) x = 748;
        if(y < 1) y = 1;
        if(y > 498) y = 498;

        gc.fillStyle = "black";
        gc.fillRect(x - 2, y - 2, 4, 4);
    }
}


