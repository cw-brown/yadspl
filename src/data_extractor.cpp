#include "data_extractor.h"

DataExtractor::DataExtractor(
        std::vector<uint8_t> * arq_buffer,
        std::vector<uint8_t> * ack_buffer,
        std::vector<SimpPacket> * data_packet_buffer) :
        ackBuffer(ack_buffer),
        arqBuffer(arq_buffer),
        dataPacketBuffer(data_packet_buffer) {

    goodPackets = new std::vector<SimpPacket>;
    goodPackets->resize(sizeof(uint8_t));

    // Set all packets in good packets to missing
    for(std::vector<SimpPacket>::iterator set = goodPackets->begin(); set != goodPackets->cend(); set++) {

        // Denote missing with control code 0b1110, which is not a valid code per the spec
        set->controlCode = 0b1110;

    }

}

void DataExtractor::processNextPacket() {

    bool matched;

    inputPos = dataPacketBuffer->begin();

    // If packet has bad checksum request repeat else acknowledge and treat as good
    if(!(inputPos->erc)) {

        arqBuffer->push_back(inputPos->sequenceNumber);
        

    } else {

        ackBuffer->push_back(inputPos->sequenceNumber);

        // put good packet to good packet buffer
        goodPackets->at(inputPos->sequenceNumber) = *inputPos;

        // We need to request repeats for missing packets
        for(std::vector<SimpPacket>::iterator check = goodPackets->begin(); check != goodPackets->cend(); check++) {

            // Missing packets marked with control code 0b1110
            if(check->controlCode == 0b1110) {

                // Check packet isn't already requested
                matched = false;
                for(std::vector<uint8_t>::iterator currArq = arqBuffer->begin(); currArq != arqBuffer->cend(); currArq++) {

                    if((*currArq) == (check->sequenceNumber)) {

                        matched = true; // Packet already arqed

                    }

                }

                // If not already in arq bufer, add it
                if(!matched) {

                    arqBuffer->push_back(std::distance(goodPackets->begin(), check));

                }

            }

        }

    }

    // Remove processed packet from buffer
    dataPacketBuffer->erase(inputPos);

}

std::vector<SimpPacket>::iterator DataExtractor::getInputPosition() {

    return inputPos;

}

uint8_t DataExtractor::getGoodPacketsContinuous() {

    uint8_t countGood;


    countGood = 0;

    for(std::vector<SimpPacket>::iterator check = goodPackets->begin(); check != goodPackets->cend(); check++) {

        if(check->controlCode == 0b1110) {

            countGood++;

        }

        // If there are ever more packets checked than good packets, then good packets must not be continuous
        if(countGood < std::distance(goodPackets->begin(), check) + 1) {

            return countGood - 1;

        }

    }

    return countGood;

}

bool DataExtractor::getGoodPacketsFull() {

    return (getGoodPacketsContinuous() == goodPackets->size());
    
}

void DataExtractor::flushGoodPacketBuffer() {

    for(std::vector<SimpPacket>::iterator toClear = goodPackets->begin(); toClear != goodPackets->cend(); toClear++) {

        toClear->controlCode = 0b1110;

    }

}