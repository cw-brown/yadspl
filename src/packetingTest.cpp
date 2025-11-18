# include <iostream>
# include <random>

# include <vector>

# include "ptransmitter.h"
# include "preceiver.h"

int main() {

    const uint8_t sAddr= 0;
    const uint8_t rAddr = 1;

    const uint8_t dLen = 2;

    const float cenFr = 1e3;

    const uint8_t mType = 0;

    const float ber = 1e-3;

    const size_t bytesToSend = 1000;

    // Send Buffers
    std::vector<bool> *  sendDataIn;
    std::vector<SimpPacket> * sendControlPs;
    std::vector<SimpPacket> * sendDataPs;
    std::vector<uint8_t> * sendData;
    std::vector<Packet> * sendPacketsToSend;

    sendData->resize(bytesToSend);

    // Receive
    std::vector<bool> *  recDataIn;
    std::vector<SimpPacket> * recControlPs;
    std::vector<SimpPacket> * recDataPs;
    std::vector<Packet> * recDataPs;
    std::vector<uint8_t> * toAck;
    std::vector<uint8_t> * toArq;
    std::vector<uint8_t> * dataReceived;
    std::vector<Packet> * recPacketsToSend;

    // Sender class
    PacketTransmitter send(sAddr, rAddr, dLen, cenFr, mType, sendDataIn, sendControlPs, sendDataPs, sendData,  sendPacketsToSend);

    // Receiver class
    PacketReceiver rec(sAddr, rAddr, dLen, cenFr, mType, recDataIn, recControlPs, recDataPs, toAck, toArq, dataReceived, recPacketsToSend);

    // Random
    std::default_random_engine generator;

    // Make random data
    std::uniform_int_distribution<uint8_t> rDat(0,7);

    // Random bit corruptions
    std::bernoulli_distribution err(ber);

    // Generate random data
    for(std::vector<uint8_t>::iterator byte = sendData->begin(); byte != sendData->end(); byte++) {

        *byte = rDat(generator);

    }

    // Transmit loop
    while(rec.getState() != REND) {

        send.tick();

        while(!sendPacketsToSend->empty()) {

            singlePacketToVector(sendPacketsToSend->front(), recDataIn);
            sendPacketsToSend->erase(sendPacketsToSend->begin());

        }

        while(!recPacketsToSend->empty()) {

            singlePacketToVector(recPacketsToSend->front(), sendDataIn);
            recPacketsToSend->erase(recPacketsToSend->begin());

        }

        rec.tick();

    }

    for(int i = 0; i < dataReceived->size(); i++) {

        std::cout << (dataReceived->at(i) == sendData->at(i));

    }

}