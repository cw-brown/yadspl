#include "preceiver.h"

PacketReceiver::PacketReceiver(
    uint8_t sender_address,
    uint8_t receiver_address,
    uint8_t data_length,
    float center_freq,
    uint8_t modulation_type,
    std::vector<bool> * in_data_buffer,
    std::vector<SimpPacket> * control_packet_buffer,
    std::vector<SimpPacket> * data_packet_buffer,
    std::vector<uint8_t> * arq_buffer,
    std::vector<uint8_t> * ack_buffer,
    std::vector<uint8_t> * receive_data,
    std::vector<Packet> * to_send) :
    senderAddress(sender_address),
    receiverAddress(receiver_address),
    centerFreq(center_freq),
    modulationType(modulation_type),
    dataLength(data_length),
    pFormer(sender_address, receiver_address),
    pReformer(in_data_buffer, control_packet_buffer, data_packet_buffer),
    dExtractor(arq_buffer, ack_buffer, data_packet_buffer),
    receiveData(receive_data),
    toSend(to_send),
    state(START),
    sequenceNumber(0) {

    dataPackets = new std::vector<Packet>;
    dataPackets->resize(sizeof(uint8_t));

}

void PacketReceiver::tick() {

    // Always process whatever has been received, also do whatever is appropriate for state
    pReformer.formPacket(false, senderAddress, receiverAddress);

    switch (state) {

        case START:
            
            // Wait for source to intiate transaction
            if(!controlPacketBuffer->empty() && ((*(controlPacketBuffer->end()-1)).controlCode == 0b1000 || (*(controlPacketBuffer->end()-1)).controlCode == 0b1001)) {

                state = WAIT_DATA;
                controlPacketBuffer->pop_back();

            }
            break;

        case WAIT_DATA:

            if(!controlPacketBuffer->empty()) {
                
                state = CONTROL;

            } else if(!dataPacketBuffer->empty()) {

                state = UNPACKET;

            }
            break;

        case CONTROL:
            
            switch ((controlPacketBuffer->back()).controlCode) {
                
                // Busy not implemented
                /*case 0b0010:
                    
                    // Busy start
                    toSend->push_back(pFormer.formBusyAcceptPacket());
                    state = BUSY;
                    break;*/

                // Frequency change not implemented
                /*case 0b0100:

                // Center frequency change
                break;*/

                case 0b0101:

                    // Modulation Change
                    if((controlPacketBuffer->back()).erc) {

                        toSend->push_back(pFormer.formGeneralChannelAcknowledge());
                        modulationType = ((controlPacketBuffer->back()).data).front();

                    }
                    state = WAIT_DATA;
                    break;

                // Transfer not implemented
                /*case 1010:
                    
                    // Transfer
                    break;*/

                case 1011:
                case 1111:

                    // Ended or Dropped
                    state = ENDING;
                    break;

                // Unknown control packet, so do nothing
                default:
                    
                    break;

                
            }

            // Remove processed control packet
            controlPacketBuffer->pop_back();
            break;

        case UNPACKET:

            // Process next data packet
            dExtractor.processNextPacket();

            if(!ackBuffer->empty()) {

                for(std::vector<uint8_t>::iterator ack = ackBuffer->begin(); ack != ackBuffer->end(); ack++) {

                    toSend->push_back(pFormer.formAcknowledgePacket(0, *ack));

                }

                ackBuffer->clear();

            }

            if(!arqBuffer->empty()) {

                for(std::vector<uint8_t>::iterator arq = arqBuffer->begin(); arq != arqBuffer->end(); arq++) {

                    toSend->push_back(pFormer.formRepeatPacket(*arq));

                }

                arqBuffer->clear();

            }

            // Need to check if we are still waiting on data in this frame
            if(dExtractor.getGoodPacketsFull()) {

                state = OUT_BUFF;

            } else{

                state = WAIT_DATA;
            }

        case OUT_BUFF:
            
            dExtractor.extractGoodPackets(receiveData);
            state = START;
            break;

        case ENDING:
            
            dExtractor.extractGoodPackets(receiveData);
            toSend->push_back(pFormer.formTransactionEndAcknowledge());
            state = END;
            break;

        case END:
        default:

            break;
    }
}

RState PacketReceiver::getState() {

    return state;
    
}

float PacketReceiver::getCenterFrequency() {

    return centerFreq;

}

uint8_t PacketReceiver::getModulationType() {

    return modulationType;

}