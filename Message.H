// ===============================================================
// Message.H 
//
// For message passing in the MPI (parallel) implementation of the
// training program.
//
// bmajoros@duke.edu
// ===============================================================
#ifndef INCL_MESSAGE_H
#define INCL_MESSAGE_H
#include "BOOM/GSL/Vector.H"


/****************************************************************
                         enum MessageType
 ****************************************************************/
enum MessageType {
  MSG_SET_INTERVAL,  // master -> slave : "use this interval"
  MSG_EVALUATE,      // master -> slave : "evaluate the likelihood"
  MSG_LIKELIHOOD,    // slave -> master : "the likelihood is L"
  MSG_TREE_NOT_FIXED,// master -> slave : "allow branch lengths to vary"
  MSG_NEXT_HIGHER_ORDER, // master -> slave : "change model order to N"
  MSG_QUIT           // master -> slave : "terminate"
};



/****************************************************************
                         class Message
 ****************************************************************/

class Message {
  void *buffer, *firstSlot;
  GSL::Vector parms;
  int numParms, bufferSize;
  long begin, end;
  double likelihood;
  int order;

  void packType(MessageType);
  void computeBufferSize();
public:
  Message(int numParms);
  ~Message();
  void *getBuffer();
  int getBufferSize();
  void resize(int numParms);

  // For sender:
  void packInterval(long begin,long end);
  void packEvaluate(const GSL::Vector &parms);
  void packQuit();
  void packLikelihood(double likelihood);
  void packTreeNotFixed();
  void packNextHigherOrder(int order,const GSL::Vector &parms);

  // For receiver:
  MessageType unpack();
  GSL::Vector &getParms();
  long getBegin();
  long getEnd();
  double getLikelihood();
  int getOrder();
};



#endif
