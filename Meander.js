"use strict";
var Spline = /** @class */ (function () {
    function Spline(xs, ys) {
        this.xs = xs;
        this.ys = ys;
        this.ks = this.getNaturalKs(new Float64Array(this.xs.length));
    }
    Spline.prototype.getIndexEnd = function (index) {
        return this.xs[index];
    };
    Spline.prototype.getIndexEndY = function (index) {
        return this.ys[index];
    };
    Spline.prototype.getNaturalKs = function (ks) {
        var n = this.xs.length - 1;
        var A = zerosMat(n + 1, n + 2);
        for (var i = 1; i < n; i++ // rows
        ) {
            A[i][i - 1] = 1 / (this.xs[i] - this.xs[i - 1]);
            A[i][i] =
                2 *
                    (1 / (this.xs[i] - this.xs[i - 1]) + 1 / (this.xs[i + 1] - this.xs[i]));
            A[i][i + 1] = 1 / (this.xs[i + 1] - this.xs[i]);
            A[i][n + 1] =
                3 *
                    ((this.ys[i] - this.ys[i - 1]) /
                        ((this.xs[i] - this.xs[i - 1]) * (this.xs[i] - this.xs[i - 1])) +
                        (this.ys[i + 1] - this.ys[i]) /
                            ((this.xs[i + 1] - this.xs[i]) * (this.xs[i + 1] - this.xs[i])));
        }
        A[0][0] = 2 / (this.xs[1] - this.xs[0]);
        A[0][1] = 1 / (this.xs[1] - this.xs[0]);
        A[0][n + 1] =
            (3 * (this.ys[1] - this.ys[0])) /
                ((this.xs[1] - this.xs[0]) * (this.xs[1] - this.xs[0]));
        A[n][n - 1] = 1 / (this.xs[n] - this.xs[n - 1]);
        A[n][n] = 2 / (this.xs[n] - this.xs[n - 1]);
        A[n][n + 1] =
            (3 * (this.ys[n] - this.ys[n - 1])) /
                ((this.xs[n] - this.xs[n - 1]) * (this.xs[n] - this.xs[n - 1]));
        return solve(A, ks);
    };
    Spline.prototype.getIndexBefore = function (target) {
        var low = 0;
        var high = this.xs.length;
        var mid = 0;
        while (low < high) {
            mid = Math.floor((low + high) / 2);
            if (this.xs[mid] < target && mid !== low) {
                low = mid;
            }
            else if (this.xs[mid] >= target && mid !== high) {
                high = mid;
            }
            else {
                high = low;
            }
        }
        if (low === this.xs.length - 1) {
            return this.xs.length - 1;
        }
        return low + 1;
    };
    Spline.prototype.at = function(targetX, segmentIndexParam) {
        var segmentIndex = segmentIndexParam ? segmentIndexParam : this.getIndexBefore(targetX);
        var normalizedX = (targetX - this.xs[segmentIndex - 1]) / (this.xs[segmentIndex] - this.xs[segmentIndex - 1]);
    
        var deltaK = this.ks[segmentIndex - 1] * (this.xs[segmentIndex] - this.xs[segmentIndex - 1]) -
                     (this.ys[segmentIndex] - this.ys[segmentIndex - 1]);
    
        var deltaKNext = -this.ks[segmentIndex] * (this.xs[segmentIndex] - this.xs[segmentIndex - 1]) +
                         (this.ys[segmentIndex] - this.ys[segmentIndex - 1]);
    
        var interpolation = (1 - normalizedX) * this.ys[segmentIndex - 1] +
                            normalizedX * this.ys[segmentIndex] +
                            normalizedX * (1 - normalizedX) * (deltaK * (1 - normalizedX) + deltaKNext * normalizedX);
    
        return interpolation;
    };
    Spline.prototype.derivativeUsingDeltaK = function(x, segmentIndex) {
        let x_i_minus_1 = this.xs[segmentIndex - 1];
        let x_i = this.xs[segmentIndex];
        let y_i_minus_1 = this.ys[segmentIndex - 1];
        let y_i = this.ys[segmentIndex];
        let k_i_minus_1 = this.ks[segmentIndex - 1];
        let k_i = this.ks[segmentIndex];

        let Delta_x = x_i - x_i_minus_1;
        let Delta_y = y_i - y_i_minus_1;

        let t = (x - x_i_minus_1) / Delta_x;
        let deltaK = k_i_minus_1 * Delta_x - Delta_y;
        let deltaKNext = -k_i * Delta_x + Delta_y;

        let derivative = (Delta_y + (1 - 2*t) * (deltaK * (1 - t) + deltaKNext * t) + t * (1 - t) * (deltaKNext - deltaK)) / Delta_x;
        return derivative;
    }
    Spline.prototype.numberSegments = () => {
        return this.xs.length;
    }
    return Spline;
}());
function solve(A, ks) {
    var m = A.length;
    var h = 0;
    var k = 0;
    while (h < m && k <= m) {
        var i_max = 0;
        var max = -Infinity;
        for (var i = h; i < m; i++) {
            var v_1 = Math.abs(A[i][k]);
            if (v_1 > max) {
                i_max = i;
                max = v_1;
            }
        }
        if (A[i_max][k] === 0) {
            k++;
        }
        else {
            swapRows(A, h, i_max);
            for (var i = h + 1; i < m; i++) {
                var f = A[i][k] / A[h][k];
                A[i][k] = 0;
                for (var j = k + 1; j <= m; j++)
                    A[i][j] -= A[h][j] * f;
            }
            h++;
            k++;
        }
    }
    for (var i = m - 1; i >= 0; i-- // rows = columns
    ) {
        var v = 0;
        if (A[i][i]) {
            v = A[i][m] / A[i][i];
        }
        ks[i] = v;
        for (var j = i - 1; j >= 0; j-- // rows
        ) {
            A[j][m] -= A[j][i] * v;
            A[j][i] = 0;
        }
    }
    return ks;
}
function zerosMat(r, c) {
    var A = [];
    for (var i = 0; i < r; i++)
        A.push(new Float64Array(c));
    return A;
}
function swapRows(m, k, l) {
    var p = m[k];
    m[k] = m[l];
    m[l] = p;
}

function distanceFromSpline(startingX, startingY, radius, x, splineIndex, spline){
    return Math.pow(x - startingX,2) + Math.pow(spline.at(x,splineIndex) - startingY,2) - Math.pow(radius,2);
}

function distanceFromSplineDerivative(startingX, startingY, x, splineIndex, spline){
    return 2 * (x - startingX + (spline.at(x,splineIndex) - startingY) * spline.derivativeUsingDeltaK(x, splineIndex));
}

function newtonsRoot(startingX, startingY, radius, initialGuess, threshold, splineIndex, spline){
    let currentGuess = initialGuess;
    let currentVal;
    let iterations = 0;
    while(Math.abs((currentVal = distanceFromSpline(startingX, startingY, radius, currentGuess, splineIndex, spline))) > threshold && iterations < 20){ // Getting closer to zero intersection
        currentGuess = currentGuess - currentVal/distanceFromSplineDerivative(startingX, startingY, currentGuess, splineIndex, spline);
        iterations++;
    }
    return currentGuess;
}

function tricube(nPoints){
    const nPointsSpacedBetween = Array(nPoints).fill().map((_,i) => -1 + 2/(nPoints - 1) * i);
    return nPointsSpacedBetween.map(val => Math.pow(1 - Math.pow(Math.abs(val),3),3));
}

function curvature_non_normal(xArr, yArr){
    const xsArr = diff_ds(xArr);
    const ysArr = diff_ds(yArr);
    const xssArr = diff_ds2(xArr);
    const yssArr = diff_ds2(yArr);
    let epsilon = 1e-7
    const kappa = xsArr.map((xsVal, i) => Math.abs(xsVal * yssArr[i] - ysArr[i] * xssArr[i]) / (Math.pow(xsVal * xsVal + ysArr[i] * ysArr[i], 1.5) + epsilon) );
    return kappa;
}

function diff_ds(numArr){
    let dsArr = new Array(numArr.length).fill(0);

    // Central difference for the bulk of the array
    for (let i = 1; i < numArr.length - 1; i++) {
        dsArr[i] = 0.5 * (numArr[i + 1] - numArr[i - 1]);
    }

    // Forward difference for the first element
    dsArr[0] = numArr[1] - numArr[0];

    // Backward difference for the last element
    dsArr[numArr.length - 1] = numArr[numArr.length - 1] - numArr[numArr.length - 2];

    return dsArr;
}

function diff_ds2(numArr){
    let dsArr = new Array(numArr.length).fill(0); // Initialize an array of zeros with the same length as x

    // Second derivative for the bulk of the array using central difference method
    for (let i = 1; i < numArr.length - 1; i++) {
        dsArr[i] = numArr[i - 1] - 2 * numArr[i] + numArr[i + 1];
    }

    // For the first and last elements, the function uses a simple difference
    // which is not a true second derivative but provides a boundary condition.
    dsArr[0] = numArr[1] - numArr[0];
    dsArr[numArr.length - 1] = numArr[numArr.length - 1] - numArr[numArr.length - 2];

    return dsArr;
}

function cumulative_distance_non_norm(xArr,yArr){
    // Initialize the array for cumulative distances with the first element as 0
    let distanceArr = [0];

    // Calculate the cumulative distance
    for (let i = 1; i < xArr.length; i++) { //indexing correct?
        let dx = xArr[i] - xArr[i - 1];
        let dy = yArr[i] - yArr[i - 1];
        let segmentLength = Math.sqrt(dx * dx + dy * dy); //pythagorean
        distanceArr.push(distanceArr[i - 1] + segmentLength);
    }

    return distanceArr;
}

function cubic_pulse(nPoints){
    const cubicPulseArr = Array(nPoints).fill()
        .map((_,i) => Math.abs(-1 + (2 / (nPoints - 1)) * i))
        .map(val => 1 - (val * val * (3 - 2 * val)))
    return cubicPulseArr;
}

function convolve_reflect(input, weights) {
    let filterSize = weights.length;
    let reversedWeights = [...weights].reverse();
    let size1 = Math.floor(filterSize / 2);
    let size2 = filterSize - size1 - 1;
    let output = new Array(input.length).fill(0.0);
    let origin = 0;

    if (filterSize % 2 === 0) {
        origin = -1;
    }

    for (let ll = 0; ll < input.length; ll++) {
        for (let jj = -size1; jj <= size2; jj++) {
            let index = ll + jj - origin;
            let weightIndex = jj + size1; 

            // Apply mode-specific boundary handling
            let value = 0;
            if (index >= 0 && index < input.length) {
                value = input[index];
            } else {
                if (index < 0) {
                    index = Math.abs(index + 1);
                } else if (index >= input.length) {
                    index = 2 * input.length - index - 1;
                }
                value = input[Math.min(Math.max(index, 0), input.length - 1)];
            }

            if (weightIndex >= 0 && weightIndex < filterSize) {
                output[ll] += value * reversedWeights[weightIndex];
            }
        }
    }

    return output;
}

function smoothingGaussian(numArr, ir){
    const kernelArr = cubic_pulse(2 * ir + 1);
    const kernelSum = kernelArr.reduce((acc, val) => acc + val, 0);
    const normalizedKernelArr = kernelArr.map(k => k / kernelSum);
    const convolution = convolve_reflect(numArr, normalizedKernelArr);
    return convolution;
}

function ccw(A, B, C) {
    return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0]);
}

function intersect(A, B, C, D) {
    return ccw(A, C, D) !== ccw(B, C, D) && ccw(A, B, C) !== ccw(A, B, D);
}

function removeLoops(xArr, yArr) {
    // Make copies of x and y to avoid modifying the original arrays
    const xArrCopy = xArr.slice();
    const yArrCopy = yArr.slice();

    let k = 0;
    while (k < xArrCopy.length - 1) {
        let removeNode = false;
        let pBreak = -1; // Store the index where the loop was found
        for (let p = k + 2; p < xArrCopy.length - 1; p++) {
            let a = [xArrCopy[k], yArrCopy[k]];
            let b = [xArrCopy[k + 1], yArrCopy[k + 1]];
            let c = [xArrCopy[p], yArrCopy[p]];
            let d = [xArrCopy[p + 1], yArrCopy[p + 1]];
            if (intersect(a, b, c, d)) {
                removeNode = true;
                pBreak = p;
                break;
            }
        }

        if (removeNode) {
            // Splice method changes the original array to remove elements,
            // and returns the removed elements as a new array

            xArrCopy.splice(k + 1, pBreak - k);
            yArrCopy.splice(k + 1, pBreak - k);
            k = pBreak + 1;
        } else {
            k += 1;
        }
    }

    return { xArr: xArrCopy, yArr: yArrCopy };
}

function interp_spline_cubic_y(xArr, yArr, nPoints=50){
    let obj = interp_spline_cubic(xArr,yArr,nPoints)
    return obj.yArr
}

function interp_spline_cubic_x(xArr, yArr, nPoints=50){
    let obj = interp_spline_cubic(xArr,yArr,nPoints)
    return obj.xArr
}

function getInitialGuess(initialX, radius, splineIndex, spline){
    let derivative = spline.derivativeUsingDeltaK(initialX, splineIndex);
    if (Number.isNaN(derivative)) {
        return initialX + radius * 0.5;
    } else {
        return initialX +  Math.sqrt(Math.pow(radius,2)/(1 + Math.pow(derivative,2))) * (spline.getIndexEnd(splineIndex) > initialX ? 1 : -1);
    }
}

function getDistance(x1,y1,x2,y2){
    return Math.sqrt(Math.pow(x1 - x2,2) + Math.pow(y1 - y2,2));
}

function interp_spline_cubic(xArr, yArr, nPoints=50) {
    const distanceArr = cumulative_distance_non_norm(xArr, yArr);
    const spline = new Spline(xArr, yArr);
    const totalLengthAdjusted = distanceArr[distanceArr.length - 1] * 1.1;
    const distancePerSegment = totalLengthAdjusted / nPoints;
    const interpolatedXArray = [];
    const interpolatedYArray = [];
    interpolatedXArray.push(xArr[0]);
    interpolatedYArray.push(yArr[0]);

    const endX = xArr[xArr.length - 1];
    const endY = yArr[yArr.length - 1];
    let currentSplineIndex = 1;
    let root;

    for (let i = 0; i < nPoints; i++) {
        if (getDistance(interpolatedXArray[i], interpolatedYArray[i], endX, endY) <= distancePerSegment) {
            interpolatedXArray.push(endX);
            interpolatedYArray.push(endY);
            break;
        }

        if (currentSplineIndex >= xArr.length) {
            break; // Exit the loop if no more valid segments are available
        }
        //Has some strange skipping behavior at higher derivatives, hence I will use the bisection method
        let initialGuess = getInitialGuess(interpolatedXArray[i], distancePerSegment, currentSplineIndex, spline);
        while (currentSplineIndex < xArr.length && ((root = newtonsRoot(interpolatedXArray[i], interpolatedYArray[i], distancePerSegment, initialGuess, 0.0001, currentSplineIndex, spline)) > Math.max(spline.getIndexEnd(currentSplineIndex), interpolatedXArray[i]) || root <= Math.min(spline.getIndexEnd(currentSplineIndex), interpolatedXArray[i]))) {
            currentSplineIndex++;
        }

        interpolatedXArray.push(root);
        interpolatedYArray.push(spline.at(root, currentSplineIndex));
    }

    interpolatedXArray[interpolatedXArray.length - 1] = endX;
    interpolatedYArray[interpolatedYArray.length - 1] = endY;
    
    return { xArr: interpolatedXArray, yArr: interpolatedYArray };
}

function manageSegments(xArr, yArr, numberOfPoints = 50) {
    function calculateSegmentLengths(xArr, yArr) {
        let lengths = [];
        for (let i = 1; i < xArr.length; i++) {
            let dx = xArr[i] - xArr[i - 1];
            let dy = yArr[i] - yArr[i - 1];
            lengths.push(Math.sqrt(dx * dx + dy * dy));
        }
        return lengths;
    }

    function calculateAverageLength(lengths) {
        let totalLength = lengths.reduce((acc, len) => acc + len, 0);
        return totalLength / lengths.length;
    }

    function findSegmentToBisect(lengths) {
        let maxLen = -Infinity;
        let index = -1;
        for (let i = 0; i < lengths.length; i++) {
            if (lengths[i] > maxLen) {
                maxLen = lengths[i];
                index = i;
            }
        }
        return index;
    }

    function bisectSegment(xArr, yArr, index) {
        if (index < 0 || index >= xArr.length - 1) {
            return;
        }

        let newX = (xArr[index] + xArr[index + 1]) / 2;
        let newY = (yArr[index] + yArr[index + 1]) / 2;

        xArr.splice(index + 1, 0, newX);
        yArr.splice(index + 1, 0, newY);
    }

    function findSegmentToRemove(lengths) {
        let minLen = Infinity;
        let index = -1;
        for (let i = 0; i < lengths.length; i++) {
            if (lengths[i] < minLen) {
                minLen = lengths[i];
                index = i;
            }
        }
        return index;
    }

    function removeSegment(xArr, yArr, index) {
        if (index <= 0 || index >= xArr.length - 1) {
            return;
        }

        xArr.splice(index, 1);
        yArr.splice(index, 1);
    }

    let lengths = calculateSegmentLengths(xArr, yArr);
    let averageLength = calculateAverageLength(lengths);
    let maxThreshold = 3 * averageLength;
    let minThreshold = averageLength / 10
    
    while(findSegmentToRemove(lengths) !== -1 && lengths[findSegmentToRemove(lengths)] < minThreshold){
        lengths = calculateSegmentLengths(xArr, yArr);
        let segmentToRemove = findSegmentToRemove(lengths);
        if (segmentToRemove !== -1) {
            removeSegment(xArr, yArr, segmentToRemove);
        }
    }

    while (xArr.length < numberOfPoints || (findSegmentToBisect(lengths) !== -1 && lengths[findSegmentToBisect(lengths)] > maxThreshold)) {
        let segmentToBisect = findSegmentToBisect(lengths);
        if (segmentToBisect !== -1) {
            bisectSegment(xArr, yArr, segmentToBisect);
        }

        lengths = calculateSegmentLengths(xArr, yArr);
    }

    while(xArr.length > numberOfPoints){
        lengths = calculateSegmentLengths(xArr, yArr);
        let segmentToRemove = findSegmentToRemove(lengths);
        if (segmentToRemove !== -1) {
            removeSegment(xArr, yArr, segmentToRemove);
        }
    }

    return { xArr, yArr };
}


function meander(xArr, yArr, ir, tangentRatio, normalRatio, nPointsMini, curvatureNormalization = 1){
    const shapeFactor = tricube(xArr.length);
    const kappa = curvature_non_normal(xArr, yArr).map(val => val/curvatureNormalization);

    const dxArr = Array(xArr.length).fill(0);
    const dyArr = Array(yArr.length).fill(0);

    for(let k = 0; k < xArr.length; k++){
        //local angle
        const angle = k < (xArr.length - 1) ? 
            Math.atan2(yArr[k + 1] - yArr[k], xArr[k + 1] - xArr[k]) :
            Math.atan2(yArr[k] - yArr[k - 1], xArr[k] - xArr[k - 1])
        //orientation
        const cp = k < (xArr.length - 2) ? cross_product(xArr[k + 1] - xArr[k], yArr[k + 1] - yArr[k], xArr[k + 2] - xArr[k], yArr[k + 2] - yArr[k]) : 1;
        //normal vector
        const normalX = Math.cos(angle);
        const normalY = Math.sin(angle);
        //tangent vector 
        const tangentX = Math.cos(angle + Math.sign(cp) * Math.PI / 2);
        const tangentY = Math.sin(angle + Math.sign(cp) * Math.PI / 2);

        dxArr[k] += tangentRatio * kappa[k] * tangentX;
        dyArr[k] += tangentRatio * kappa[k] * tangentY;

        dxArr[k] += normalRatio * kappa[k] * normalX;
        dyArr[k] += normalRatio * kappa[k] * normalY;
    }
    // backup length before deformation
    const previousLengthArr = cumulative_distance_non_norm(xArr, yArr);
    const previousLength = previousLengthArr[previousLengthArr.length - 1];

    //apply deformation
    const smoothDxArr = smoothingGaussian(dxArr, ir);
    const smoothDyArr = smoothingGaussian(dyArr, ir);

    const deformedXArr = xArr.map((xVal,i) => xVal + smoothDxArr[i] * shapeFactor[i]);
    const deformedYArr = yArr.map((yVal,i) => yVal + smoothDyArr[i] * shapeFactor[i]);

    const newLengthArr = cumulative_distance_non_norm(deformedXArr, deformedYArr);
    const newLength = newLengthArr[previousLengthArr.length - 1];

    const interpolatedObject = manageSegments(deformedXArr,deformedYArr);

    return removeLoops(interpolatedObject.xArr, interpolatedObject.yArr);
}

function cross_product(ux, uy, vx, vy) {
    return uy * vx - ux * vy;
}

function meanderMidpoints(xIn, yIn, amplitude, itterations=1) {
    let xArr = [...xIn];
    let yArr = [...yIn];
    for (let i = 0; i < itterations; i++) {
        const xmArr = []
        const ymArr = []
        const cp = cross_product(xArr[1] - xArr[0], yArr[1] - yArr[0], xArr[2] - xArr[0], yArr[2] - yArr[0])
        let cpSign = Math.sign(cp)
        for (let k = 0; k < xArr.length - 1; k++) {
            const angle = Math.atan2(yArr[k + 1] - yArr[k], xArr[k + 1] - xArr[k]);
            const dist = Math.hypot(xArr[k + 1] - xArr[k], yArr[k + 1] - yArr[k]);

            let xmid = 0.5 * (xArr[k + 1] + xArr[k])
            let ymid = 0.5 * (yArr[k + 1] + yArr[k])

            const da = cpSign * Math.PI / 2
            xmid += amplitude * dist * Math.cos(angle + da)
            ymid += amplitude * dist * Math.sin(angle + da)
            xmArr.push(xArr[k]);
            xmArr.push(xmid);
            ymArr.push(yArr[k]);
            ymArr.push(ymid);

            cpSign *= -1;
        }

        xmArr.push(xArr[xArr.length - 1])
        ymArr.push(yArr[yArr.length - 1])

        xArr = [...xmArr];
        yArr = [...ymArr];
    }
    return {xArr: xArr, yArr: yArr};
}

function initializeStateObject() {
   // seed = 6
   // npt = 10
   const numberOfInterpolatedPoints = 100;

    // --- first start with a set of random nodes

    // rng = np.random.default_rng(seed)
    const rngArr = Array(3).fill().map(e => Math.random() * .05 + 0.45);

    const xArr = [0,0.25,0.5,0.75,1];// make this stink less
    const yArr = [0.5,...rngArr,0.5];

    const stateObject = meanderMidpoints(xArr,yArr,0.1,2)
    // --- try to make a continuous path with these nodes using a nearest
    // --- neighbor search

    // x, y = tools.nearest_neighbor_search_2d(x, y)
    return interp_spline_cubic(stateObject.xArr, stateObject.yArr, numberOfInterpolatedPoints);
}

window.onload = () => {
    const canvas = document.getElementById("canvas");
    const ctx = canvas.getContext("2d");
    const stateObject = initializeStateObject();
    const metaObject = initializeMetaObject();
    window.sim = metaObject;
    window.state = stateObject;

    resizeCanvas(ctx, stateObject);
    window.addEventListener('resize', () => resizeCanvas(ctx, stateObject));
    animate(ctx, stateObject);
    scheduleNextUpdate(stateObject, metaObject);

    document.getElementById('playButton').addEventListener('click', () => {
        if (!metaObject.isUpdating) {
            metaObject.isUpdating = true;
            scheduleNextUpdate(stateObject, metaObject);
        }
    });

    document.getElementById('pauseButton').addEventListener('click', () => {
        metaObject.isUpdating = false;
    });
}

function setCanvasFullScreen(canvas) {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
}

function resizeCanvas(ctx, stateObject) {
    setCanvasFullScreen(ctx.canvas)
    drawState(ctx, stateObject);
}

function drawState(ctx, stateObject) { // need to fix this to not convert, also last point missing
    ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height); 

    drawStateOriginal(ctx,stateObject)
    drawStateOriginal2(ctx,stateObject)

    ctx.strokeStyle = '#21b0fe'; // Set color to red
    ctx.lineWidth = 10;
    ctx.beginPath();

    // Ensure there are at least two points to draw
    if (stateObject.xArr.length < 2) return;

    // Convert stateObject coordinates to canvas coordinates
    const points = stateObject.xArr.map((xVal, i) => ({
        x: xVal * ctx.canvas.width,
        y: stateObject.yArr[i] * ctx.canvas.height
    }));

    // Move to the first point
    ctx.moveTo(points[0].x, points[0].y);

    // Loop through the points, using the midpoint technique for quadratic curves
    for (var i = 1; i < points.length - 2; i++) {
        var xc = (points[i].x + points[i + 1].x) / 2;
        var yc = (points[i].y + points[i + 1].y) / 2;
        ctx.quadraticCurveTo(points[i].x, points[i].y, xc, yc);
    }

    // Curve through the last two points
    // This check ensures that we don't attempt this if there are less than 3 points
    if (points.length > 2) {
        ctx.quadraticCurveTo(
            points[i].x, points[i].y,
            points[i + 1].x, points[i + 1].y
        );
    }
    ctx.stroke();
}

function updateState(stateObject, onUpdateComplete) {
    const newStateObject = meander(
        stateObject.xArr,
        stateObject.yArr,
        5,
        0.02,
        0.01,
        50,
        (1 / 0.002))
    stateObject.xArr = newStateObject.xArr;
    stateObject.yArr = newStateObject.yArr;
    onUpdateComplete();
}

function initializeMetaObject() {
    return {isUpdating:true, updateInterval:10, nextUpdateDelay:10};
}

function scheduleNextUpdate(stateObject, metaObject) {
    setTimeout(() => {
        if (!metaObject.isUpdating) return; //Stop updates
        const updateStartTime = Date.now();
        updateState(stateObject, () => {
            const updateEndTime = Date.now();
            const elapsedTime = updateEndTime - updateStartTime;
            metaObject.nextUpdateDelay = Math.max(0, metaObject.updateInterval - elapsedTime);
            scheduleNextUpdate(stateObject, metaObject);
        });
    }, metaObject.nextUpdateDelay);
}

function animate(ctx, stateObject) {
    drawState(ctx, stateObject);
    requestAnimationFrame(() => animate(ctx, stateObject));
}

function drawStateOriginal(ctx, stateObject) {
    ctx.strokeStyle = '#fed700'; // White color
    ctx.lineWidth = 10;
    ctx.beginPath();
    stateObject.xArr.forEach((xVal, i) => {
        var canvasX = xVal * ctx.canvas.width;
        var canvasY = stateObject.yArr[i] * ctx.canvas.height + 50;
        if (i == 0) ctx.moveTo(canvasX, canvasY);
        else ctx.lineTo(canvasX, canvasY); // Calculate control point and switch to bezier or quadratic
    });
    ctx.stroke();
    ctx.closePath();
}

function drawStateOriginal2(ctx, stateObject) {
    ctx.strokeStyle = '#fe218b'; // White color
    ctx.lineWidth = 10;
    ctx.beginPath();
    stateObject.xArr.forEach((xVal, i) => {
        var canvasX = xVal * ctx.canvas.width;
        var canvasY = stateObject.yArr[i] * ctx.canvas.height -50;
        if (i == 0) ctx.moveTo(canvasX, canvasY);
        else ctx.lineTo(canvasX, canvasY); // Calculate control point and switch to bezier or quadratic
    });
    ctx.stroke();
    ctx.closePath();
}
