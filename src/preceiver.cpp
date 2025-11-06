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

            if(!dataPacketBuffer->empty()) {

                state = UNPACKET;
            }
            break;

        case UNPACKET:

            dExtractor.processNextPacket();
            while(!ackBuffer->empty()) {
            }
            while(!arqBuffer->empty()){
            }
            break;
        
        default:

            break;
    }
}