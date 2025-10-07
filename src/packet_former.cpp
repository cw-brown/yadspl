#include "packer_former.h"

PacketFormer::PacketFormer(
        std::ring<uint8_t> * input_buffer,
        std::ring<Packet> * output_buffer,
        uint32_t sender_address,
        uint32_t reciever_address,
        uint8_t data_length) :
        inputBuffer(input_buffer),
        outputBuffer(output_buffer),
        senderAddress(sender_address),
        recieverAddress(reciever_address),
        dataLength(data_length) {

    // Iterators for input and output buffers
    this->inputPosition = input_buffer->begin();
    this->outputPosition = output_buffer->begin();

}

void PacketFormer::formNextPacket() {

    static Packet tempPacket;

    // Reset packet contents
    (tempPacket.data).clear();
    
    // Data packet type
    tempPacket.controlCode = 0b0000;

    // Data length
    tempPacket.dataLength = dataLength;

    // Add sender and reciever addresses
    tempPacket.senderAddress = senderAddress;
    tempPacket.recieverAddress = recieverAddress;

    // Add current packet number and increment it
    tempPacket.sequenceNum = sequenceNum;
    sequenceNum++;

    // Form data section

}