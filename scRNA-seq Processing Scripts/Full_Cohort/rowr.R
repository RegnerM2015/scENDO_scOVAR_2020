#' Row-Based Functions for R Objects
#'
#' Provides utilities which interact with all R objects as
#' if they were arranged in rows.  It allows more consistent and predictable 
#' output to common functions, and generalizes a number of utility functions to
#' to be failsafe with any number and type of input objects.
#' @name rowr
#' @docType package
#' @importFrom methods as
NULL

#' Vectorize a scalar function to work on any R object.
#' 
#' Robust alternative to \code{\link{Vectorize}} function that accepts any function with two 
#' or more arguments.  Returns a function that will work an arbitrary number of vectors, lists or 
#' data frames, though output may be unpredicatable in unusual applications.  The 
#' results are also intended to be more intuitive than Vectorize.
#' 
#' @param fun a two or more argument function
#' @param type like \code{MARGIN} in \code{\link{apply}}, except that \code{c(1,2)} is
#'   represented as a \code{3} instead.  By default, will \code{Reduce} single dimensional
#'   data handle everything else row-wise.
#' @export
#' @examples
#' vectorize(`+`)(c(1,2,3))
#' vectorize(sum)(c(1,2,3),c(1,2,3))
#' # Compare these results to Vectorize, which does not vectorize sum at all.
#' Vectorize(sum)(c(1,2,3),c(1,2,3))
#' # Across data frame columns.
#' df<-data.frame(a=c(1,2,3),b=c(1,2,3))
#' vectorize(sum)(df$a,df$b)
#' # Once again, Vectorize gives a different result
#' Vectorize(sum)(df$a,df$b)
#' # Any combination of vectors, lists, matrices, or data frames an be used.
#' vectorize(`+`)(c(1,2,3),list(1,2,3),cbind(c(1,2,3)))
vectorize<-function(fun,type=NULL)
{
  function(...)
  {
    cols<-cbind.fill(...)
    if(is.null(type))
      if(dim(cols)[2]<2)
        type=2
      else
        type=1
      if(type==3)
        margin=c(1,2)
      else
        margin=type
      if(type %in% c(1,2,3))
        apply(cols,margin,function (x) Reduce(fun,unlist(x)))
      else
        Reduce(fun,unlist(cols))
  }  
}


#'Pads an object to a desired length, either with replicates of itself or another repeated object.
#'
#'@param x an R object
#'@param length.out the desired length of the final output
#'@param fill R object to fill empty rows in columns below the max size.  If unspecified, repeats input rows in the same way as \code{cbind}.
#'@param preserveClass determines whether to return an object of the same class as the original argument.  Otherwise, returns a matrix.
#'@export
#'@examples
#'buffer(c(1,2,3),20)
#'buffer(matrix(c(1,2,3,4),nrow=2),20)
#'buffer(list(1,2,3),20)
#'df<-data.frame(as.factor(c('Hello','Goodbye')),c(1,2))
#'buffer(df,5)
#'buffer((factor(x=c('Hello'))),5)
buffer<-function(x,length.out=len(x),fill=NULL,preserveClass=TRUE)
{
  xclass<-class(x)
  input<-lapply(vert(x),unlist)
  results<-as.data.frame(lapply(input,rep,length.out=length.out))
  if(length.out>len(x) && !is.null(fill))
  {
    results<-t(results)
    results[(length(unlist(x))+1):length(unlist(results))]<-fill
    results<-t(results)
  }
  if(preserveClass)
    results<-as2(results,xclass)
  return(results)   
}


#' Combine arbitrary data types, filling in missing rows.
#' 
#' Robust alternative to \code{\link{cbind}} that fills missing values and works
#' on arbitrary data types.  Combines any number of R objects into a single matrix, with each input
#' corresponding to the greater of 1 or ncol.  \code{cbind} has counterintuitive
#' results when working with lists, cannot handle certain inputs of differing
#' length, and does not allow the fill to be specified.
#' 
#' @param ... any number of R data objects
#' @param fill R object to fill empty rows in columns below the max size.  If unspecified, repeats input rows in the same way as \code{cbind}. Passed to \code{\link{buffer}}.
#' @export
#' @examples
#' cbind.fill(c(1,2,3),list(1,2,3),cbind(c(1,2,3)))
#' cbind.fill(rbind(1:2),rbind(3:4))
#'df<-data.frame(a=c(1,2,3),b=c(1,2,3))
#' cbind.fill(c(1,2,3),list(1,2,3),cbind(c('a','b')),'a',df)
#' cbind.fill(a=c(1,2,3),list(1,2,3),cbind(c('a','b')),'a',df,fill=NA)
cbind.fill<-function(...,fill=NULL)
{
  inputs<-list(...)
  inputs<-lapply(inputs,vert)
  maxlength<-max(unlist(lapply(inputs,len)))
  bufferedInputs<-lapply(inputs,buffer,length.out=maxlength,fill,preserveClass=FALSE)
  return(Reduce(cbind.data.frame,bufferedInputs))
}

#'Allows row indexing without knowledge of dimensionality or class.
#'
#'@param data any \code{R} object
#'@param rownums indices of target rows
#'@export
#'@examples
#'rows(c('A','B','C'),c(1,3))
#'rows(list('A','B','C'),c(1,3))
#'df<-data.frame(a=c(1,2,3),b=c(1,2,3))
#'rows(df,3)
rows <- function(data,rownums)
{
  #result<-data[rownums]
  if(is.null(dim(data)))
  {
    result<-data[rownums]
  }
  else
  {
    result<-data[rownums,]
  }
  #result<-ifelse(is.null(dim(data)),data[c(rownums)],data[c(rownums),])
  return((result))
}

#'Allows finding the 'length' without knowledge of dimensionality.
#'
#'@param data any \code{R} object
#'@export
#'@examples
#'len(list(1,2,3))
#'len(c(1,2,3,4))
#'df<-data.frame(a=c(1,2,3),b=c(1,2,3))
#'len(df)
len <- function(data)
{
  result<-ifelse(is.null(nrow(data)),length(data),nrow(data))
  return(result)
}


#' A more versatile form of the T-SQL \code{coalesce()} function.
#' 
#' Little more than a wrapper for \code{\link{vectorize}}, allows for 
#' duplication of SQL coalesce functionality, certain types of if-else 
#' statements, and \code{\link{apply}}/\code{\link{Reduce}} combinations.
#' 
#' @param ... an arbitrary number of \code{R} objects
#' @param fun a two argument function that returns an atomic value
#' @export
#' @examples
#' coalesce(c(NA,1,2))
#' coalesce(c(NA,1,2),c(3,4,NA))
#' df<-data.frame(a=c(NA,2,3),b=c(1,2,NA))
#' coalesce(df$a,df$b)
#' # Or even just:
#' coalesce(df)
#' # Coalesce can actually use any comparison.  For example, instead of non-NA
#' # values it could find the max in each row:
#' cbind(EuStockMarkets,Max=coalesce(EuStockMarkets,fun=function (x,y) if (x>y) x else y))
coalesce<-function(...,fun=(function (x,y) if(!is.na(x)) x else y))
{
  
  FUN=match.fun(fun)
  vectorize(FUN)(...)
}

#'A more versatile form of the T-SQL \code{count()} function.
#'
#'Implementation of T-SQL \code{count} and Excel \code{COUNTIF} functions. 
#'Shows the total number of elements in any number of data objects altogether or
#'that match a condition.
#'
#'@param ... an arbitrary number of \code{R} objects
#'@param condition a 1 argument condition
#'@export
#'@examples
#'count(c(NA,1,2))
#'count(c(NA,1,2),is.na)
#'count(c(NA,1,2),list('A',4),cbind(1,2,3))
#'count(c(NA,1,2),list('A',4),cbind(1,2,3),condition=is.character)
count<-function(...,condition=(function (x) TRUE))
{
  data<-c(...)
  result<-sum(sapply(data, function (x) if(condition(x)) 1 else 0))
  return(result)
}

#'Applies a function over a rolling window on any data object.
#'
#'Simple generalized alternative to \code{\link[zoo]{rollapply}} in package
#'\code{\link[zoo]{zoo}} with the advantage that it works on any type of data
#'structure (vector, list, matrix, etc) instead of requiring a \code{zoo}
#'object.
#'
#'@param data any \code{R} object
#'@param fun the function to evaluate
#'@param window window width defining the size of the subset available to the
#'  fun at any given point
#'@param minimum minimum width of the window.  Will not return results if the
#'  window is truncated below this value at the end of the data set
#'@param align whether to align the window right or left
#'@param ... additional arguments to pass to \code{fun}
#'@export
#'@examples
#'rollApply(1:100,sum,minimum=2,window=2)
#'rollApply(c(1,2,3),sum)
#'##6 5 3
#'rollApply(c(1,2,3,4,5,6,7,8,9),sum)
#'##45 44 42 39 35 30 24 17  9
#'rollApply(c(1,2,3,4,5,6,7,8,9),sum,window=2)
#'##3  5  7  9 11 13 15 17  9
#'rollApply(list(1,2,3,4,5,6,7,8,9),function(x) sum(unlist(x)),window=2,minimum=2)
#'##3  5  7  9 11 13 15 17
#'cbind(women,Rolling3=rollApply(women,fun=function(x) mean(x$weight),window=3,align='right'))
#'
rollApply <- function(data,fun,window=len(data),minimum=1,align='left',...)
{
  if(minimum>len(data))
    return()
  FUN=match.fun(fun)
  if (align=='left')
    result<-sapply(1:(len(data)-minimum+1),function (x) FUN(rows(data,x:(min(len(data),(x+window-1)))),...))
  if (align=='right')
    result<-sapply(minimum:len(data),function (x) FUN(rows(data,max(1,x-window+1):x),...))
  return(result)
}

#'Applies a function row-wise on any data object.
#'
#'Essentially functions as a \code{MARGIN=1} \code{\link{apply}} apply but also
#'works on data objects without 2 dimensions such as lists and vectors.
#'
#'@param data any \code{R} object
#'@param fun the function to evaluate
#'@param ... additional arguments to pass to \code{fun}
#'@export
#'@examples
#'rowApply(list(1,2,3),function (x) sum(unlist(x)))
#'df<-data.frame(a=c(1,2,3),b=c(1,2,3))
#'rowApply(df,sum)
rowApply<-function(data,fun,...)
{
  sapply(1:len(data),function (x) fun(rows(data,x),...))
}


#'Inserts a matrix into another matrix.
#'
#'Inserts a matrix or data frame into another matrix or data frame.  The new
#'rows are placed together at the row index specified.
#'@param existing table to insert into
#'@param insert rows to insert
#'@param r index at which to insert
#'@export
#'@examples
#'df1<-data.frame(a=c(1,2,3),b=c(1,2,3),c=c(1,2,3))
#'insertRows(df1,data.frame(list('a','a','a')),5)
#'insertRows(df1,data.frame(list('a','a','a')),4)
#'insertRows(df1,data.frame(list('a','a','a')),3)
#'insertRows(df1,data.frame(list('a','a','a')),2)
#'insertRows(df1,data.frame(list('a','a','a')),1)
#'insertRows(df1,df1,3)
insertRows <- function(existing, insert, r) {
  colnames(insert)<-colnames(existing)
  result<-rbind(existing,insert)
  sizeA<-len(existing)
  sizeB<-len(insert)
  order1<-c(seq_len(sizeA),rep.int(r,sizeB))
  order2<-c(rep.int(sizeB+1,sizeA),seq_len(sizeB))
  return(result[order(order1,order2),,drop=FALSE])
}

#'A more robust form of the R \code{\link{as}} function.
#'
#' Alternative to \code{as} that allows any data object to be converted to any other.  
#'
#'@param object any \code{R} object
#'@param class the name of the class to which \code{object} should be coerced
as2<-function(object,class)
{
  object<-as.matrix(object)
  if(class=='factor')
    return(as.factor(as.character(object)))
  if(class=='data.frame')
    return(as.data.frame(object))
  else
    return(as(object,class))
}

vert<-function(object)
{
  #result<-as.data.frame(cbind(as.matrix(object)))
  if(is.list(object))
    object<-cbind(object)
  object<-data.frame(object)
  
  return(object)
}